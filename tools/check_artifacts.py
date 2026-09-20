"""Check portable documentation links, manuscript tables, and public artifacts."""

from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import argparse
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]


class Document(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.tables, self.row, self.cell = [], [], None, None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.tables.append([])
        if tag == "tr":
            self.row = []
        if tag in {"th", "td"}:
            self.cell = ""
        self.links.extend(v for k, v in attrs if k in {"href", "src"} and v)

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data

    def handle_endtag(self, tag):
        if tag in {"th", "td"}:
            self.row.append(self.cell)
            self.cell = None
        if tag == "tr":
            self.tables[-1].append(self.row)
            self.row = None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, default=ROOT / "site")
    args = parser.parse_args()
    site = args.site.resolve()
    checked = 0
    for path in site.rglob("*.html"):
        doc = Document()
        doc.feed(path.read_text())
        for raw in doc.links:
            parts = urlsplit(raw)
            if parts.scheme or parts.netloc or not parts.path:
                continue
            # An explicitly selected externally hosted MathJax installation is
            # outside this generated directory; no paper/supplement link is.
            if parts.path.startswith("/") and parts.path.endswith("tex-svg.js"):
                continue
            dest = (
                (site / parts.path.lstrip("/"))
                if parts.path.startswith("/")
                else path.parent / unquote(parts.path)
            )
            assert dest.is_file() or dest.is_dir(), (path.name, raw)
            checked += 1
    source = (ROOT / "paper/manuscript.md").read_text()
    doc = Document()
    doc.feed((site / "paper/manuscript.html").read_text())
    blocks = re.findall(
        r"^\|[^\n]+\|\n\|(?:[ \t]*:?-+:?[ \t]*\|)+[ \t]*\n(?:\|[^\n]*\|[ \t]*(?:\n|$))+",
        source,
        re.M,
    )
    normalize = lambda x: re.sub(r"\s+", " ", x).strip()
    assert len(blocks) == len(doc.tables) == 28
    for block, actual in zip(blocks, doc.tables):
        lines = block.strip().splitlines()
        expected = [
            [normalize(c) for c in line.strip().strip("|").split("|")]
            for line in [lines[0]] + lines[2:]
        ]
        assert expected == [[normalize(c) for c in row] for row in actual], expected[0]
    # Distribution must have no source references to a research workstation,
    # private mail exports, or a Tailscale host. Synthetic CSV identifiers remain.
    for folder in ["src", "docs", "paper", "results", "examples", "tests", "tools"]:
        for path in (ROOT / folder).rglob("*"):
            if (
                path.suffix not in {".py", ".md", ".json", ".csv", ".cjs"}
                or path.name == Path(__file__).name
            ):
                continue
            text = path.read_text()
            assert not re.search(
                r"/Users/|git_astra|svn_workingcopy|svn_phdthesis|100\.74\.231|Mail Export|draft_methods\.pdf",
                text,
            ), path
    for item in json.loads((ROOT / "docs/paper_provenance.json").read_text()):
        assert (
            hashlib.sha256((ROOT / item["artifact"]).read_bytes()).hexdigest()
            == item["sha256"]
        ), item["artifact"]
    print(
        json.dumps(
            {
                "passed": True,
                "local_links": checked,
                "manuscript_tables": len(doc.tables),
                "all_table_cells_match": True,
                "portable_source_scan": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
