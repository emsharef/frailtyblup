"""Compile the canonical LaTeX paper with Tectonic."""

from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
SOURCE = PAPER / "manuscript.tex"
OUTPUT = PAPER / "paper.pdf"


def main() -> None:
    tectonic = shutil.which("tectonic")
    if tectonic is None:
        raise SystemExit(
            "Tectonic is required: https://tectonic-typesetting.github.io/"
        )

    with tempfile.TemporaryDirectory(prefix="frailtyblup-paper-") as temporary:
        output_dir = Path(temporary)
        subprocess.run(
            [tectonic, "--outdir", str(output_dir), SOURCE.name],
            cwd=PAPER,
            env={**os.environ, "SOURCE_DATE_EPOCH": "1789862400", "TZ": "UTC"},
            check=True,
        )
        compiled = output_dir / "manuscript.pdf"
        if not compiled.is_file() or compiled.stat().st_size == 0:
            raise SystemExit("Tectonic completed without producing manuscript.pdf")
        temporary_pdf = OUTPUT.with_suffix(".pdf.tmp")
        shutil.copyfile(compiled, temporary_pdf)
        temporary_pdf.replace(OUTPUT)

    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
