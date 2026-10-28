# Method overview

The browser calculator accepts **one** H5AD at a time, with raw unnormalized counts in an explicitly selected layer and gene symbols matched to a chosen bundled fixed reference. Cells can optionally be filtered by an H5AD observation column.

The browser loads the frozen scientific Python engine through Pyodide. Gene expression is scaled without mean centering and per-cell fixed-reference non-negative matrix factorization reconstructs program weights. The usage output normalizes nonzero-cell weights to sum to 100%. The calculator does **not** relearn spectra from uploaded data.

The three references are human all-cell (18 programs, 2,999 genes), human myeloid (14, 2,278 genes) and human-derived mouse ortholog-mapped myeloid (14, 1,897 mouse genes). The references and scientific worker are protected by SHA-256 checks.

Downloads contain usage_percent.tsv, raw_w.tsv, run_qc.json and provenance.json. cell_qc.tsv is optional through Advanced options; QC runs even if the file is omitted.

## Limitations
- This release supports one H5AD file per calculation.
- Raw-count validation is not yet comprehensive: a misleadingly named layer may contain normalized counts.
- Active input cap: 512 MiB. Large 10–15 GiB experimental runs did not establish an approved memory-safe cap.
- Research use, not clinical interpretation.

The underlying software and reference licenses remain in LICENSE, site/python/LICENSE and site/assets/references/LICENSE.
