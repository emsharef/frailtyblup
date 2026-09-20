"""Build a portable static documentation site without reading external research files."""

import argparse
import html
from pathlib import Path
import re
import shutil
from urllib.parse import urlsplit, urlunsplit
import markdown

ROOT = Path(__file__).resolve().parents[1]


def render(source, out, mathjax):
    text = source.read_text()
    formulas = []

    def preserve(match):
        display = match.group(0).startswith("$$")
        tex = match.group(1) if display else match.group(2)
        token = f"MATHPLACEHOLDER{len(formulas):05d}TOKEN"
        formulas.append((token, tex, display))
        return "\n\n" + token + "\n\n" if display else token

    protected = re.sub(r"\$\$(.*?)\$\$|\\\((.*?)\\\)", preserve, text, flags=re.S)
    body = markdown.markdown(protected, extensions=["tables", "fenced_code", "toc"])
    for token, tex, display in formulas:
        if display:
            body = body.replace(
                "<p>" + token + "</p>",
                '<div class="math-display">$$' + html.escape(tex) + "$$</div>",
            )
        else:
            body = body.replace(
                token, '<span class="math-inline">\\(' + html.escape(tex) + "\\)</span>"
            )

    def local_link(match):
        prefix, url, suffix = match.groups()
        parts = urlsplit(html.unescape(url))
        if not parts.scheme and not parts.netloc and parts.path.endswith(".md"):
            path = parts.path[:-3] + ".html"
            if path.endswith("README.html"):
                path = path[: -len("README.html")] + "index.html"
            url = urlunsplit(
                (parts.scheme, parts.netloc, path, parts.query, parts.fragment)
            )
        return prefix + html.escape(url, quote=True) + suffix

    body = re.sub(r'((?:href|src)=")([^"]+)(")', local_link, body)
    body = body.replace("<table>", '<div class="table-wrap"><table>').replace(
        "</table>", "</table></div>"
    )
    if source.name == "web-version.md":
        body = body.replace(
            "<table>\n<thead>\n<tr>\n<th>Analysis cohort</th>",
            '<table class="cohort-table">\n<thead>\n<tr>\n<th>Analysis cohort</th>',
        )
        sections = [
            (ident, label)
            for ident, label in re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
            if label != "Abstract"
        ]
        if sections:
            contents = (
                '<section class="manuscript-contents"><h2>Contents</h2><ul>'
                + "".join(
                    f'<li><a href="#{ident}">{label}</a></li>'
                    for ident, label in sections
                )
                + "</ul></section>"
            )
            body = body.replace(
                '<h2 id="' + sections[0][0] + '">',
                contents + '<h2 id="' + sections[0][0] + '">',
                1,
            )
        body = re.sub(
            r'<p><strong>Table (?:6c|6d|6e|10c|10d)\..*?</p>\s*<div class="table-wrap"><table>.*?</table></div>',
            lambda m: '<div class="keep-table">' + m.group(0) + "</div>",
            body,
            flags=re.S,
        )
    rel = source.relative_to(ROOT)
    target = out / (
        "index.html" if rel == Path("README.md") else str(rel.with_suffix(".html"))
    )
    prefix = "../" * len(rel.parent.parts)
    heading = re.search(r"^# (.+)$", text, re.M)
    title = heading.group(1) if heading else source.stem
    nav = " · ".join(
        f'<a href="{prefix}{url}">{label}</a>'
        for url, label in [
            ("index.html", "FrailtyBLUP"),
            ("paper/paper.pdf", "Paper"),
            ("docs/data.html", "Data"),
            ("docs/api.html", "API"),
            ("docs/examples.html", "Examples"),
            ("docs/reproducibility.html", "Reproducibility"),
        ]
    )
    math_config = r"window.MathJax={tex:{inlineMath:[['\\(','\\)']],displayMath:[['$$','$$']]},svg:{fontCache:'local'},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};"
    klass = (
        "review-manuscript thesis-extension-manuscript"
        if source.name == "web-version.md"
        else "package-docs"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><link rel="stylesheet" href="{prefix}assets/style.css"><script>{math_config}</script><script defer src="{html.escape(mathjax,quote=True)}"></script></head><body class="{klass}"><nav>{nav}</nav><main>{body}</main></body></html>'
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "site")
    parser.add_argument(
        "--mathjax-url",
        default="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js",
    )
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "assets").mkdir(exist_ok=True)
    shutil.copyfile(ROOT / "docs/style.css", out / "assets/style.css")
    sources = (
        list(ROOT.glob("*.md"))
        + list((ROOT / "docs").glob("*.md"))
        + list((ROOT / "paper").glob("*.md"))
    )
    for source in sources:
        render(source, out, args.mathjax_url)
    for folder in ["paper", "results", "examples", "src", "tools", "docs"]:
        for source in (ROOT / folder).rglob("*"):
            if (
                not source.is_file()
                or source.suffix in {".md", ".pyc"}
                or any(x in source.parts for x in ["__pycache__", "output"])
            ):
                continue
            target = out / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
    for name in ["LICENSE", "CITATION.cff", "pyproject.toml"]:
        shutil.copyfile(ROOT / name, out / name)
    print(f"Built {len(sources)} HTML pages in {out}")


if __name__ == "__main__":
    main()
