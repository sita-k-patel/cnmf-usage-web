# cNMF Usage Calculator

Browser-based calculator for fixed-reference single-cell program usage in glioma datasets.

**Website:** https://cnmf-usage.tymillerlab.org
**Source:** https://github.com/sita-k-patel/cnmf-usage-web

## What it does

Upload one AnnData `.h5ad` containing raw unnormalized counts in a named layer, select a fixed reference and optional cell populations, then calculate per-cell program usage directly in the browser. The calculator does not learn new gene programs.

### Bundled references

| Reference | Programs | Genes |
|---|---:|---:|
| Human glioma, all cells | 18 | 2,999 |
| Human glioma, myeloid cells | 14 | 2,278 |
| Mouse glioma, myeloid cells (human-derived) | 14 | 1,897 |

The mouse reference is the 14 human glioma myeloid programs mapped to mouse orthologous genes, **not** de novo mouse-derived programs.

## Input and results

- Accepts **one H5AD** per run.
- Requires raw, unnormalized, nonnegative expression counts in a named layer (not merely `.X`).
- Gene names must match reference gene symbols; select a suitable gene-symbol source when necessary.
- Optional cell selection uses observation metadata.
- Current browser input limit: **512 MiB**. Larger test runs were experimental and are not an advertised capacity.
- The raw-count validator is not yet comprehensive. The layer name alone cannot prove the contents are raw counts.

The default ZIP contains `usage_percent.tsv`, `raw_w.tsv`, `run_qc.json`, and `provenance.json`. Under Advanced options, enable `cell_qc.tsv` for per-cell fit/QC information; QC calculations always run.

This tool is for research, not clinical diagnosis or treatment.

## Reference and software citation

Human myeloid programs: Miller TE, El Farran CA, Couturier CP, et al. *Programs, origins and immunomodulatory functions of myeloid cells in glioma.* Nature 640, 1072–1082 (2025). https://doi.org/10.1038/s41586-025-08633-8

Mouse ortholog-mapped spectra derive from those human programs; an associated cross-species manuscript is under review. Do not cite it as published.

Scientific package: [cnmf-usage](https://github.com/cwru-non-academic/cnmf-usage).
Software license: MIT (Ty Miller Lab), see `LICENSE` and `site/python/LICENSE`. The distributed reference spectra have separate attribution/licensing requirements; see `site/assets/references/LICENSE`.

## Development

Run locally:

```bash
cd site
python3 -m http.server 8770 --bind 127.0.0.1
```

Open http://127.0.0.1:8770/.

Run frozen-source integrity and lightweight UI checks from the repository root:

```bash
python3 tools/check_integrity.py
python3 -m unittest discover -s tests/ui -p 'test_*.py' -v
```

These checks do **not** replace full browser numerical parity tests. Before publishing science-affecting changes, run the private three-fixture D1/B6.3 browser regression and review the outputs.

## Publishing

Cloudflare Pages should publish **only `site/`**, not the repository root. The production branch is `main`. Keep private datasets, results ZIPs, and internal handoffs outside this public GitHub repository.
See [methods](docs/METHODS.md), [validation](docs/VALIDATION.md), and [release procedure](docs/RELEASE.md).
