"""Make generated PDF link annotations independent of a local preview server."""

import argparse
from pathlib import Path
import posixpath
from urllib.parse import urlsplit, urlunsplit
import pymupdf


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--base-url", required=True)
    args = parser.parse_args()
    doc = pymupdf.open(args.pdf)
    changed = 0
    for page in doc:
        for link in page.get_links():
            uri = link.get("uri", "")
            if uri.startswith(args.base_url.rstrip("/") + "/"):
                parts = urlsplit(uri[len(args.base_url.rstrip("/")) + 1 :])
                relative = posixpath.relpath(parts.path, "paper")
                link["uri"] = urlunsplit(
                    ("", "", relative, parts.query, parts.fragment)
                )
                page.update_link(link)
                changed += 1
    doc.set_metadata(
        {
            "title": "Linear Prediction and Frailty Inference for Clustered Recurrent Events",
            "subject": "Manuscript 0.12: public companion repository",
            "creator": "FrailtyBLUP documentation",
        }
    )
    content = doc.tobytes(garbage=4, deflate=True)
    doc.close()
    args.pdf.write_bytes(content)
    check = pymupdf.open(args.pdf)
    assert all(
        args.base_url not in link.get("uri", "")
        for p in check
        for link in p.get_links()
    )
    print(f"Made {changed} PDF links portable; {len(check)} pages.")


if __name__ == "__main__":
    main()
