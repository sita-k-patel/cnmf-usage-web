# Deployment and release

1. Confirm the `main` production site currently works. Back up Cloudflare Pages configuration.
2. Review changed files and exclude all patient data, private handoffs, temporary ZIPs and credentials.
3. Run `python3 tools/check_integrity.py` and `python3 -m unittest discover -s tests/ui -p 'test_*.py' -v`.
4. Run the real D1/B6.3 H5AD browser regression, both ZIP variants and negative input tests.
5. For the one-time migration, **change Cloudflare Pages build output directory from `.` to `site` only when the `site/` copy is already deployed on `main`**. Avoid removing root copies until the new production site works.
6. Confirm `https://cnmf-usage.tymillerlab.org` loads, receives JavaScript and Pyodide assets, runs a permitted synthetic H5AD, and downloads matching results.
7. Only then remove duplicate root website paths and obsolete development pin files in a separate PR. Monitor deployment and keep a rollback commit.

Do not include a site version JSON in this release.
