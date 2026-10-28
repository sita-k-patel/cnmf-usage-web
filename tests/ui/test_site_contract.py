"""Lightweight source-level UI contract tests. These do not replace real-browser science regression."""
import unittest
from pathlib import Path
from html.parser import HTMLParser

ROOT = Path(__file__).resolve().parents[2]
class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = []
        self.links = []
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a: self.ids.append(a["id"])
        if tag == "a" and "href" in a: self.links.append(a["href"])

class TestUI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.html = (ROOT / "site/index.html").read_text()
        cls.css = (ROOT / "site/css/site.css").read_text()
        cls.page = Page()
        cls.page.feed(cls.html)
    def test_unique_ids(self):
        self.assertEqual(len(self.page.ids), len(set(self.page.ids)))
    def test_core_controls(self):
        for item in ("file-input","layer-select","gene-select","include-cell-qc",
                     "calc","cite-open","cite-open-start","cite","cite-copy"):
            self.assertIn(item,self.page.ids)
    def test_reference_metadata(self):
        for phrase in ("2,999 genes","2,278 genes","1,897 genes",
                       "human-derived","mouse orthologous genes"):
            self.assertIn(phrase,self.html)
    def test_navigation_and_css(self):
        self.assertIn("https://github.com/sita-k-patel/cnmf-usage-web",self.page.links)
        self.assertIn("Citation",self.html)
        self.assertIn("#282D37",self.css)
        self.assertIn("#86AEFF",self.css)
    def test_bundle_paths(self):
        for rel in ("site/js/ui.js","site/js/backend-worker.js",
                    "site/js/worker.mjs","site/python/engine.py",
                    "site/assets/references/LICENSE"):
            self.assertTrue((ROOT / rel).is_file(),rel)
if __name__ == "__main__":
    unittest.main()
