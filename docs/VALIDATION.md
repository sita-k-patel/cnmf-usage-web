# Validation status

## Evidence so far

The frozen scientific engine and worker match the checksums in `integrity/frozen.json`. Prior D1 development recorded three small fixtures passing with nine TSVs byte-identical to B6.3 browser baselines, **before** UI adapter changes.

A subsequent 2,000-cell public-UI calculation produced 18 program columns; its ZIP passed structural and numerical consistency checks. This is not a substitute for byte-for-byte baseline comparison.

## Remaining verification

- Rerun all three small D1/B6.3 fixtures against the deployed UI adapter and compare original scientific TSV values.
- Test the 4-file default and 5-file optional-QC ZIP variants.
- Run a synthetic known-answer fixture with an independent expected output.
- Test wrong species/genes, unsupported layers and normalized or negative matrices.
- Record Chrome console and download failures; only publish verified claims.

The current source-level UI and file-integrity workflows cannot establish scientific correctness. Experimental 10–15 GiB calculations completed but did not pass the memory-acceptance criteria, so the advertised cap remains 512 MiB.

## Local browser regression — 2026-10-28

Synthetic fixture: 02_gene_symbols_in_var_names.h5ad
Reference: human_adult_glioma_allcell_k18
Input: 160 cells, 2,999 shared genes, counts layer

Independent comparison: browser Pyodide engine versus
cnmf-usage Python CLI on the same synthetic input.

- NMF converged in both implementations at 20 iterations.
- Cell IDs, ordering, reference SHA-256 and program columns matched.
- raw_w.tsv: 0/2880 mismatches; maximum absolute difference 3.64e-12.
- usage_percent.tsv: 0/2880 mismatches; maximum absolute difference 9.95e-14.
- Optional cell_qc.tsv: see separate QC comparison results.
- Browser checks passed for multiple gene columns, optional QC,
  missing symbols/layer, negative counts, CSC rejection and zero-count cells.
