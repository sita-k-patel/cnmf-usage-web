"""Two-pass, chunked back-calculation engine for cnmf-usage.

ONE engine for every place it runs: the Terra tests, the local app and the
browser (Pyodide). It does NOT import scanpy, anndata or cnmf_usage.

Locked numerical contract (validated bit-identical to cnmf_usage 0.1.0):
  shared genes  -> exact match, input gene order
  gene order    -> indices sorted in each row (sum_duplicates), as scanpy does
  Pass 1        -> per-gene sum and sum of squares over ALL cells; SD (n - 1)
  Pass 2        -> x / SD (scaling rule from the scaling probe), X @ H.T
  solver        -> sklearn _update_cdnmf_fast, W = 0 start, ONE global stop rule
  usage         -> W / row_sum * 100
  QC            -> same formulas as cnmf_usage.core

Options that do NOT change the math: chunk size, cache, cell filter (masks),
cell-ID prefixes, virtual copies (benchmarks only), progress callbacks,
output sinks.
"""

from __future__ import annotations

import hashlib
import io
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import non_negative_factorization
from sklearn.decomposition._cdnmf_fast import _update_cdnmf_fast
from sklearn.utils import check_array

ENGINE_VERSION = "0.2.1"

# --------------------------------------------------------------------------
# Small helpers copied from cnmf_usage.core (identical logic)
# --------------------------------------------------------------------------


def ordered_sha256(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


class StreamingSha256:
    """Same result as ordered_sha256, but fed in chunks."""

    def __init__(self):
        self.d = hashlib.sha256()

    def update(self, values):
        for value in values:
            encoded = str(value).encode("utf-8")
            self.d.update(len(encoded).to_bytes(8, "big"))
            self.d.update(encoded)

    def hexdigest(self):
        return self.d.hexdigest()


def fit_error_summary(values):
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    keys = ["min", "p25", "p50", "p75", "p95", "max"]
    if finite.size == 0:
        return {k: None for k in keys}
    pct = np.percentile(finite, [0, 25, 50, 75, 95, 100])
    return {k: float(v) for k, v in zip(keys, pct)}


def relative_fit_error(x_norm_sq, x_h, w, h):
    """Same as cnmf_usage.core._relative_fit_error, from streamed pieces."""
    h_h = h @ h.T
    model_norm_sq = np.einsum("ij,jk,ik->i", w, h_h, w, optimize=True)
    cross_term = np.einsum("ij,ij->i", w, x_h, optimize=True)
    residual_sq = x_norm_sq + model_norm_sq - 2.0 * cross_term
    residual_sq = np.maximum(residual_sq, 0.0)
    denominator = np.sqrt(x_norm_sq)
    result = np.full(x_norm_sq.shape[0], np.nan, dtype=np.float64)
    valid = denominator > 0
    result[valid] = np.sqrt(residual_sq[valid]) / denominator[valid]
    return result


def library_versions():
    import platform
    import scipy
    import sklearn

    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
            "pandas": pd.__version__, "h5py": getattr(h5py, "__version__", None),
            "engine": ENGINE_VERSION}


# --------------------------------------------------------------------------
# Reference
# --------------------------------------------------------------------------


@dataclass
class Reference:
    id: str
    spectra: pd.DataFrame  # genes x programs
    meta: dict = field(default_factory=dict)

    @classmethod
    def from_export(cls, ref_dir, ref_id):
        """Exported from the package environment (Terra tests)."""
        ref_dir = Path(ref_dir)
        z = np.load(ref_dir / f"{ref_id}.npz", allow_pickle=False)
        meta = json.loads((ref_dir / f"{ref_id}.json").read_text())
        spectra = pd.DataFrame(
            z["values"],
            index=pd.Index([str(g) for g in z["genes"]]),
            columns=[str(p) for p in z["programs"]],
        )
        return cls(id=ref_id, spectra=spectra, meta=meta)

    @classmethod
    def from_tsv_bytes(cls, raw: bytes, meta: dict):
        """Same loading as cnmf_usage.references.load_reference (with SHA-256 check)."""
        actual = hashlib.sha256(raw).hexdigest()
        if actual != meta["sha256"]:
            raise RuntimeError(
                f"Reference '{meta['id']}' failed its SHA-256 check. "
                f"Expected {meta['sha256']}, found {actual}. The analysis stopped.")
        spectra = pd.read_table(io.BytesIO(raw), sep="\t", header=0, index_col=0)
        return cls(id=meta["id"], spectra=spectra, meta=meta)


# --------------------------------------------------------------------------
# Chunked reader for .h5ad count layers (h5py only, no anndata)
# --------------------------------------------------------------------------


def _attr_str(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def _read_strings(ds):
    """Read index strings from legacy arrays or AnnData nullable-string groups.

    Cell and gene IDs must not contain missing values. No counts are changed.
    """
    if isinstance(ds, h5py.Group):
        enc = _attr_str(ds.attrs.get("encoding-type", ""))
        if enc == "nullable-string-array":
            values = _read_strings(ds["values"])
            mask = np.asarray(ds["mask"][()], dtype=bool)
            if mask.shape != values.shape:
                raise ValueError(f"Malformed nullable strings at {ds.name}: mask shape differs")
            if np.any(mask):
                raise ValueError(f"Missing cell/gene IDs are not allowed at {ds.name}")
            return values
        if enc == "categorical":
            cats = _read_strings(ds["categories"])
            codes = np.asarray(ds["codes"][()])
            if np.any(codes < 0) or np.any(codes >= len(cats)):
                raise ValueError(f"Missing or invalid categorical cell/gene IDs at {ds.name}")
            return cats[codes]
        raise ValueError(f"Unsupported string index encoding '{enc}' at {ds.name}")
    if h5py.check_string_dtype(ds.dtype) is not None:
        return np.asarray(ds.asstr()[()], dtype=object)
    return np.asarray([_attr_str(v) for v in ds[()]], dtype=object)


def _read_index(group):
    name = _attr_str(group.attrs["_index"]) if "_index" in group.attrs else "_index"
    return _read_strings(group[name])


def read_column(group, name):
    """One obs/var column as strings, including nullable and categorical arrays."""
    node = group[name]
    if isinstance(node, h5py.Group):
        enc = _attr_str(node.attrs.get("encoding-type", ""))
        if enc == "nullable-string-array":
            values = _read_strings(node["values"])
            mask = np.asarray(node["mask"][()], dtype=bool)
            if mask.shape != values.shape:
                raise ValueError(f"Malformed nullable strings at {node.name}: mask shape differs")
            values[mask] = "nan"  # match old metadata-column missing-value handling
            return values
        if enc != "categorical":
            raise ValueError(f"Unsupported column encoding '{enc}' at {node.name}")
        cats = _read_strings(node["categories"])
        codes = np.asarray(node["codes"][()])
        if np.any(codes >= len(cats)):
            raise ValueError(f"Invalid categorical codes at {node.name}")
        out = np.empty(codes.shape[0], dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = "nan"
        return out
    cat_attr = node.attrs.get("categories", None)  # anndata < 0.8
    if cat_attr is not None:
        cats = _read_strings(group.file[cat_attr])
        codes = node[()]
        return np.asarray([cats[c] if c >= 0 else "nan" for c in codes], dtype=object)
    if h5py.check_string_dtype(node.dtype) is not None:
        return _read_strings(node)
    return np.asarray([str(v) for v in node[()]], dtype=object)


def column_names(group):
    if "column-order" in group.attrs:
        order = [_attr_str(c) for c in group.attrs["column-order"]]
        if order:
            return order
    idx = _attr_str(group.attrs.get("_index", "_index"))
    return [k for k in group.keys() if k not in (idx, "__categories")]


class H5Counts:
    """Reads rows of one count layer from an .h5ad file, chunk by chunk.

    `source` is a path or a binary file-like object (the browser passes one).
    """

    def __init__(self, source, layer, gene_col=None, name=None):
        self.path = name or str(source)
        self.f = h5py.File(source, "r")
        node = self.f["X"] if layer == "X" else self.f["layers"][layer]
        self.node = node
        if isinstance(node, h5py.Group):
            enc = _attr_str(node.attrs.get("encoding-type", ""))
            self.shape = tuple(int(x) for x in node.attrs["shape"])
            self.dtype = node["data"].dtype
            if enc == "csr_matrix":
                self.fmt = "csr"
                self.indptr = node["indptr"][()].astype(np.int64)
                self.data_ds = node["data"]
                self.ind_ds = node["indices"]
            elif enc == "csc_matrix":
                self.fmt = "csc"
                try:
                    from anndata.io import sparse_dataset
                except ImportError:
                    try:
                        from anndata.experimental import sparse_dataset
                    except ImportError as e:  # browser: no anndata
                        raise ValueError(
                            "The counts layer is stored by gene (CSC). This version "
                            "reads only cell-ordered (CSR) layers.") from e
                self.sds = sparse_dataset(node)
            else:
                raise ValueError(f"Unsupported sparse encoding: {enc}")
        else:
            self.fmt = "dense"
            self.shape = tuple(node.shape)
            self.dtype = node.dtype
        self.obs_names = _read_index(self.f["obs"])
        if gene_col:
            self.var_names = read_column(self.f["var"], gene_col)
        else:
            self.var_names = _read_index(self.f["var"])
        self.bytes_read = 0

    @property
    def n_obs(self):
        return self.shape[0]

    def rows(self, a, b):
        """Rows a..b-1 as a CSR matrix in the stored dtype (all genes)."""
        if self.fmt == "csr":
            s, e = int(self.indptr[a]), int(self.indptr[b])
            data = self.data_ds[s:e]
            ind = self.ind_ds[s:e]
            ptr = self.indptr[a : b + 1] - s
            self.bytes_read += data.nbytes + ind.nbytes
            return sparse.csr_matrix((data, ind, ptr), shape=(b - a, self.shape[1]))
        if self.fmt == "csc":
            m = sparse.csr_matrix(self.sds[a:b])
            self.bytes_read += m.data.nbytes + m.indices.nbytes
            return m
        arr = self.node[a:b, :]
        self.bytes_read += arr.nbytes
        return sparse.csr_matrix(arr)

    def close(self):
        self.f.close()


class MultiInput:
    """Several files (optional cell filter, optional virtual copies) read as ONE population."""

    def __init__(self, sources, layer, reps=1, masks=None, id_prefixes=None,
                 gene_col=None, names=None):
        names = names or [None] * len(sources)
        self.readers = [H5Counts(s, layer, gene_col=gene_col, name=nm)
                        for s, nm in zip(sources, names)]
        self.reps = int(reps)
        self.masks = masks or [None] * len(self.readers)
        self.prefixes = id_prefixes or [""] * len(self.readers)
        for r, m in zip(self.readers, self.masks):
            if m is not None and len(m) != r.n_obs:
                raise ValueError(f"{r.path}: cell filter length does not match the file")
        per = [int(np.count_nonzero(m)) if m is not None else r.n_obs
               for r, m in zip(self.readers, self.masks)]
        self.n_obs = sum(per) * self.reps
        self.dtype = self.readers[0].dtype
        for r in self.readers[1:]:
            if r.dtype != self.dtype:
                raise ValueError("All input files must use the same count dtype.")

    def iter_chunks(self, chunk_size):
        """Yield (global_start, reader_index, rep, a, b, sel); sel = selected local rows or None."""
        g = 0
        for rep in range(self.reps):
            for ri, r in enumerate(self.readers):
                m = self.masks[ri]
                for a in range(0, r.n_obs, chunk_size):
                    b = min(a + chunk_size, r.n_obs)
                    if m is None:
                        yield g, ri, rep, a, b, None
                        g += b - a
                    else:
                        sel = np.nonzero(m[a:b])[0]
                        if len(sel):
                            yield g, ri, rep, a, b, sel
                            g += len(sel)

    def n_chunks(self, chunk_size):
        return sum(1 for _ in self.iter_chunks(chunk_size))

    def rows(self, ri, a, b, sel):
        m = self.readers[ri].rows(a, b)
        return m if sel is None else m[sel]

    def cell_ids(self, ri, rep, a, b, sel=None):
        ids = self.readers[ri].obs_names[a:b]
        if sel is not None:
            ids = ids[sel]
        p = self.prefixes[ri]
        if p:
            ids = np.asarray([f"{p}{x}" for x in ids], dtype=object)
        if self.reps == 1:
            return ids
        return np.asarray([f"{x}-r{rep}" for x in ids], dtype=object)

    def base_ids(self):
        out = []
        for ri, r in enumerate(self.readers):
            ids = r.obs_names
            if self.masks[ri] is not None:
                ids = ids[self.masks[ri]]
            p = self.prefixes[ri]
            out.append(np.asarray([f"{p}{x}" for x in ids], dtype=object) if p else ids)
        return np.concatenate(out)

    @property
    def bytes_read(self):
        return sum(r.bytes_read for r in self.readers)

    def close(self):
        for r in self.readers:
            r.close()


# --------------------------------------------------------------------------
# Scaling (same math as scanpy.pp.scale(zero_center=False))
# --------------------------------------------------------------------------

# var_mode:   "sq_in"  square in the stored dtype | "sq_f64" square in float64
# scale_mode: "div_f64" float64(x)/SD | "div_in" round(float64(x)/SD) to stored dtype
#             "mul_recip_f64" | "mul_recip_in"
# Validated for cnmf-usage 0.1.0 (scanpy 1.11.5) with float32 counts: sq_in + div_in.
VAR_MODES = ("sq_in", "sq_f64")
SCALE_MODES = ("div_f64", "div_in", "mul_recip_f64", "mul_recip_in")
VALIDATED_SCALING = {"var_mode": "sq_in", "scale_mode": "div_in"}


def chunk_sums(sub, var_mode):
    G = sub.shape[1]
    data = sub.data
    s = np.bincount(sub.indices, weights=data.astype(np.float64), minlength=G)
    if var_mode == "sq_in" and data.dtype.kind == "f":
        sq = (data * data).astype(np.float64)
    else:
        d64 = data.astype(np.float64)
        sq = d64 * d64
    q = np.bincount(sub.indices, weights=sq, minlength=G)
    return s, q


def std_from_sums(s, q, n):
    mean = s / n
    var = q / n - mean**2
    var *= n / (n - 1)
    std = np.sqrt(var)
    std[std == 0] = 1
    return std


def scale_chunk(sub, std, scale_mode):
    data = sub.data
    sd = std.take(sub.indices, mode="clip")
    in_float = data.dtype.kind == "f"
    if scale_mode == "div_f64":
        out = data.astype(np.float64) / sd
    elif scale_mode == "div_in":
        out = data.astype(np.float64) / sd
        if in_float:
            out = out.astype(data.dtype).astype(np.float64)
    elif scale_mode == "mul_recip_f64":
        out = data.astype(np.float64) * (1 / std).take(sub.indices, mode="clip")
    elif scale_mode == "mul_recip_in":
        r = (1 / std).take(sub.indices, mode="clip")
        if in_float:
            tmp = data.copy()
            tmp *= r
            out = tmp.astype(np.float64)
        else:
            out = data.astype(np.float64) * r
    else:
        raise ValueError(scale_mode)
    return sparse.csr_matrix((out, sub.indices.copy(), sub.indptr.copy()), shape=sub.shape)


# --------------------------------------------------------------------------
# Solver
# --------------------------------------------------------------------------


def solve_fixed_h(XHt, HHt, tol, max_iter, on_round=None):
    """Same loop as sklearn _fit_coordinate_descent with update_H=False."""
    n, K = XHt.shape
    W = np.zeros((n, K), dtype=np.float64)
    perm = np.arange(K, dtype=np.intp)
    n_iter = 0
    violation_init = None
    for n_iter in range(1, max_iter + 1):
        violation = 0.0
        violation += _update_cdnmf_fast(W, HHt, XHt, perm)
        if on_round is not None:
            on_round(n_iter)
        if n_iter == 1:
            violation_init = violation
        if violation_init == 0:
            break
        if violation / violation_init <= tol:
            break
    return W, n_iter


def build_h(reference, shared_genes):
    h = reference.spectra.loc[shared_genes].T.to_numpy(dtype=np.float64)
    Ht = check_array(h.T, order="C")
    HHt = np.dot(Ht.T, Ht)
    return h, Ht, HHt


def shared_genes_for(var_names, reference):
    ref_genes = set(reference.spectra.index)
    return [str(g) for g in var_names if str(g) in ref_genes]


# --------------------------------------------------------------------------
# Main two-pass run
# --------------------------------------------------------------------------


def _noop(*_a, **_k):
    pass


TSV_KW = dict(sep="\t", index=True, index_label="cell_id", na_rep="NaN", float_format="%.17g")


def run_twopass(
    inputs, layer, reference, out_dir=None, *,
    chunk_size=1000, var_mode="sq_in", scale_mode="div_in", cache=False, reps=1,
    tol=1e-4, max_iter=1000, nmf_parameters=None, save_arrays=True, write_tsv=True,
    masks=None, id_prefixes=None, gene_col=None, names=None,
    outputs=("usage", "raw_w", "cell_qc"), sink=None,
    stage=_noop, progress=_noop,
):
    """sink(name, text): receives each output file chunk by chunk (browser).
    Without a sink, files are written to out_dir."""
    out_dir = Path(out_dir) if out_dir is not None else None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
    info = {}

    stage("open")
    mi = MultiInput(inputs, layer, reps=reps, masks=masks, id_prefixes=id_prefixes,
                    gene_col=gene_col, names=names)
    first = mi.readers[0]
    shared = shared_genes_for(first.var_names, reference)
    if not shared:
        raise ValueError("No shared genes between the input and the reference.")
    positions = []
    for r in mi.readers:
        vn = pd.Index([str(v) for v in r.var_names])
        if vn.has_duplicates:
            raise ValueError(f"{r.path}: duplicate gene IDs")
        pos = vn.get_indexer(shared)
        if (pos < 0).any():
            raise ValueError(f"{r.path}: some shared genes are missing in this file")
        positions.append(pos)
    n = mi.n_obs
    if n < 2:
        raise ValueError("At least two cells are required.")
    G = len(shared)
    h, Ht, HHt = build_h(reference, shared)
    K = h.shape[0]
    n_chunks = mi.n_chunks(chunk_size)

    # ---------------- Pass 1 ----------------
    stage("pass1")
    s_tot = np.zeros(G)
    q_tot = np.zeros(G)
    genes_detected = np.empty(n, dtype=np.int64)
    counts_sum = None
    cache_list = [] if cache else None
    nnz_all = nnz_shared = 0
    id_hash = StreamingSha256()
    if pd.Index(mi.base_ids()).has_duplicates:
        raise ValueError("Input contains duplicate cell IDs (use file-name prefixes).")
    for k, (g0, ri, rep, a, b, sel) in enumerate(mi.iter_chunks(chunk_size)):
        full = mi.rows(ri, a, b, sel)
        m = full.shape[0]
        nnz_all += full.nnz
        if full.nnz and (not np.isfinite(full.data).all()):
            raise ValueError("Counts contain NaN or Inf.")
        if full.nnz and (full.data < 0).any():
            raise ValueError("Counts contain negative values.")
        sub = full[:, positions[ri]]
        del full
        sub.sum_duplicates()  # same storage order as the package (locked)
        nnz_shared += sub.nnz
        s, q = chunk_sums(sub, var_mode)
        s_tot += s
        q_tot += q
        genes_detected[g0 : g0 + m] = np.asarray((sub > 0).sum(axis=1)).ravel().astype(np.int64)
        cs = np.asarray(sub.sum(axis=1)).ravel()
        if counts_sum is None:
            counts_sum = np.empty(n, dtype=cs.dtype)
        counts_sum[g0 : g0 + m] = cs
        id_hash.update(mi.cell_ids(ri, rep, a, b, sel))
        if cache:
            cache_list.append(sub)
        progress("pass1", k + 1, n_chunks)
    std = std_from_sums(s_tot, q_tot, n)

    # ---------------- Pass 2 ----------------
    stage("pass2")
    XHt = np.empty((n, K), dtype=np.float64)
    x_norm_sq = np.empty(n, dtype=np.float64)
    if cache:
        g0 = 0
        for k, sub in enumerate(cache_list):
            xs = scale_chunk(sub, std, scale_mode)
            m = sub.shape[0]
            XHt[g0 : g0 + m] = xs @ Ht
            x_norm_sq[g0 : g0 + m] = np.asarray(xs.multiply(xs).sum(axis=1)).ravel()
            g0 += m
            progress("pass2", k + 1, n_chunks)
        cache_list = None
    else:
        for k, (g0, ri, rep, a, b, sel) in enumerate(mi.iter_chunks(chunk_size)):
            sub = mi.rows(ri, a, b, sel)[:, positions[ri]]
            sub.sum_duplicates()
            m = sub.shape[0]
            xs = scale_chunk(sub, std, scale_mode)
            XHt[g0 : g0 + m] = xs @ Ht
            x_norm_sq[g0 : g0 + m] = np.asarray(xs.multiply(xs).sum(axis=1)).ravel()
            progress("pass2", k + 1, n_chunks)

    # ---------------- Fit ----------------
    stage("solve")
    t0 = time.perf_counter()
    W, n_iter = solve_fixed_h(XHt, HHt, tol, max_iter,
                              on_round=lambda i: progress("solve", i, max_iter))
    info["solve_seconds"] = time.perf_counter() - t0
    info["seconds_per_round"] = info["solve_seconds"] / max(n_iter, 1)

    stage("qc")
    raw_w_sum = W.sum(axis=1)
    fit_err = relative_fit_error(x_norm_sq, XHt, W, h)
    del XHt, x_norm_sq
    zero_usage = raw_w_sum == 0

    # ---------------- Outputs ----------------
    stage("write")
    programs = list(reference.spectra.columns)
    out_bytes = {}
    if write_tsv:
        out_bytes = _write_outputs_chunked(
            mi, out_dir, sink, outputs, W, programs, genes_detected, counts_sum,
            raw_w_sum, fit_err, zero_usage, chunk_size, progress, n_chunks)
    if save_arrays and out_dir is not None:
        ids_all = np.concatenate(
            [mi.cell_ids(ri, rep, a, b, sel) for _, ri, rep, a, b, sel in mi.iter_chunks(chunk_size)])
        with np.errstate(divide="ignore", invalid="ignore"):
            usage = W / W.sum(axis=1, keepdims=True) * 100.0
        np.savez(out_dir / "arrays.npz", cell_ids=ids_all.astype(str), raw_w=W, usage=usage,
                 reference_genes_detected=genes_detected, reference_counts_sum=counts_sum,
                 raw_w_sum=raw_w_sum, relative_fit_error=fit_err,
                 zero_usage_cell=zero_usage, std=std)

    ref_gene_set = set(shared)
    spectra = reference.spectra
    missing = [str(g) for g in spectra.index if str(g) not in ref_gene_set]
    mass = spectra.loc[shared].sum(axis=0) / spectra.sum(axis=0)
    run_qc = {
        "engine": "twopass",
        "engine_version": ENGINE_VERSION,
        "library_versions": library_versions(),
        "reference_id": reference.id,
        "reference_sha256": reference.meta.get("sha256"),
        "reference_version": reference.meta.get("reference_version"),
        "counts_layer": layer,
        "input_files": [r.path for r in mi.readers],
        "input_dtype": str(mi.dtype),
        "n_cells": int(n),
        "n_input_genes": int(len(first.var_names)),
        "n_reference_genes": int(spectra.shape[0]),
        "n_shared_genes": G,
        "gene_coverage_fraction": G / spectra.shape[0],
        "missing_reference_genes": missing,
        "program_mass_retained": {p: float(mass[p]) for p in spectra.columns},
        "nmf_n_iter": int(n_iter),
        "nmf_max_iter": int(max_iter),
        "nmf_converged": bool(n_iter < max_iter),
        "zero_usage_cells": int(zero_usage.sum()),
        "relative_fit_error": fit_error_summary(fit_err),
        "ordered_cell_ids_sha256": id_hash.hexdigest(),
        "ordered_shared_genes_sha256": ordered_sha256(shared),
        "nmf_parameters": nmf_parameters,
        "scaling": {"var_mode": var_mode, "scale_mode": scale_mode},
        "chunk_size": int(chunk_size),
        "n_input_files": len(inputs),
        "cell_filter_used": any(m is not None for m in mi.masks),
        "cell_id_prefixes": list(mi.prefixes) if any(mi.prefixes) else None,
        "virtual_copies": int(reps),
        "cache_shared_genes": bool(cache),
    }
    qc_text = json.dumps(run_qc, indent=2)
    if sink is not None:
        sink("run_qc.json", qc_text)
    if out_dir is not None:
        (out_dir / "run_qc.json").write_text(qc_text)
    info.update(
        n_cells=int(n), n_programs=int(K), n_shared_genes=G, nmf_n_iter=int(n_iter),
        nnz_all_genes=int(nnz_all), nnz_shared_genes=int(nnz_shared),
        shared_nnz_fraction=float(nnz_shared / nnz_all) if nnz_all else None,
        bytes_read_uncompressed=int(mi.bytes_read), input_dtype=str(mi.dtype),
        output_bytes=out_bytes,
        cell_id_chars_mean=float(np.mean([len(str(x)) for x in first.obs_names[:10000]])),
    )
    mi.close()
    return info


def _write_outputs_chunked(mi, out_dir, sink, outputs, W, programs, det, csum, wsum, ferr,
                           zero, chunk_size, progress, n_chunks):
    """Same TSV format as cnmf_usage.cli._write_table, written chunk by chunk."""
    names = {"usage": "usage_percent.tsv", "raw_w": "raw_w.tsv", "cell_qc": "cell_qc.tsv"}
    handles = {}
    if sink is None:
        handles = {k: open(out_dir / names[k], "w", newline="") for k in outputs}
    nbytes = {names[k]: 0 for k in outputs}

    def emit(key, df, header):
        buf = io.StringIO()
        df.to_csv(buf, header=header, **TSV_KW)
        text = buf.getvalue()
        nbytes[names[key]] += len(text.encode("utf-8"))
        if sink is not None:
            sink(names[key], text)
        else:
            handles[key].write(text)

    first = True
    for k, (g0, ri, rep, a, b, sel) in enumerate(mi.iter_chunks(chunk_size)):
        ids = mi.cell_ids(ri, rep, a, b, sel)
        sl = slice(g0, g0 + len(ids))
        idx = pd.Index(ids, name="cell_id")
        w = W[sl]
        if "usage" in outputs:
            with np.errstate(divide="ignore", invalid="ignore"):
                u = w / w.sum(axis=1, keepdims=True) * 100.0
            emit("usage", pd.DataFrame(u, index=idx, columns=programs), first)
        if "raw_w" in outputs:
            emit("raw_w", pd.DataFrame(w, index=idx, columns=programs), first)
        if "cell_qc" in outputs:
            emit("cell_qc", pd.DataFrame({
                "reference_genes_detected": det[sl], "reference_counts_sum": csum[sl],
                "raw_w_sum": wsum[sl], "relative_fit_error": ferr[sl],
                "zero_usage_cell": zero[sl]}, index=idx), first)
        first = False
        progress("write", k + 1, n_chunks)
    for fh in handles.values():
        fh.close()
    return nbytes


# --------------------------------------------------------------------------
# Simple chunking (comparison only: each chunk scaled and fitted alone)
# --------------------------------------------------------------------------


def run_naive_chunks(inputs, layer, reference, out_dir, *, chunk_size, var_mode,
                     scale_mode, nmf_parameters, stage=_noop):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stage("open")
    mi = MultiInput(inputs, layer)
    shared = shared_genes_for(mi.readers[0].var_names, reference)
    positions = [pd.Index([str(v) for v in r.var_names]).get_indexer(shared)
                 for r in mi.readers]
    h, _, _ = build_h(reference, shared)
    params = dict(nmf_parameters)
    params.pop("n_components", None)
    params.pop("random_state", None)
    stage("fit_chunks")
    Ws, ids, iters = [], [], []
    for g0, ri, rep, a, b, sel in mi.iter_chunks(chunk_size):
        sub = mi.rows(ri, a, b, sel)[:, positions[ri]]
        s, q = chunk_sums(sub, var_mode)
        std = std_from_sums(s, q, sub.shape[0])
        xs = scale_chunk(sub, std, scale_mode)
        w, _, it = non_negative_factorization(
            xs, W=None, H=h, n_components=h.shape[0], random_state=None, **params)
        Ws.append(w)
        iters.append(int(it))
        ids.append(mi.cell_ids(ri, rep, a, b, sel))
    W = np.vstack(Ws)
    with np.errstate(divide="ignore", invalid="ignore"):
        usage = W / W.sum(axis=1, keepdims=True) * 100.0
    np.savez(out_dir / "arrays.npz", cell_ids=np.concatenate(ids).astype(str),
             raw_w=W, usage=usage)
    (out_dir / "run_qc.json").write_text(json.dumps(
        {"engine": "naive_chunks", "chunk_size": chunk_size,
         "nmf_n_iter_per_chunk": iters, "n_cells": int(W.shape[0])}, indent=2))
    mi.close()
    return {"n_cells": int(W.shape[0]), "nmf_n_iter_per_chunk": iters}
