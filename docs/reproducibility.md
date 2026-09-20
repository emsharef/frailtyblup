# Reproducibility and the paper

## What this repository reproduces directly

Installation, all synthetic examples, the validated public interval interface, the constrained estimator, numerical regression tests, and documentation builds require no private files. `tests/reference.json` contains expected scores, derivatives and predictions evaluated by the original manuscript implementation on synthetic inputs. `docs/kernel_provenance.json` records the source identities. The default packaged coupled search uses the same equations, scaling and candidate-selection rule.

The paper's canonical [LaTeX source](../paper/manuscript.tex) includes all proofs, simulation tables and application results; the compiled [PDF](../paper/paper.pdf) is checked in for readers. Supporting simulation/coverage summaries, simulated per-history interval diagnostics, aggregate clinical coefficient/Jacobian summaries, and figures are in `results/`. The [artifact manifest](paper_provenance.json) identifies original and public-copy hashes. The public LaTeX source changes links and explains this repository's reproduction scope; mathematical displays and source table rows from manuscript 0.12 are preserved.

## Historical study versus standalone examples

The paper's broad study has 4,200 attempted scenario fits on 3,800 distinct histories. It uses application-informed covariate and follow-up profiles. Reconstructing those original event histories requires the underlying clinical design inputs. Those individual records are not distributed here. The new examples use an entirely synthetic generator; they should not be described as reproductions of the 4,200-fit study or evidence for new statistical claims.

The archived scenario manifest preserves cases, rates, seeds and run settings. Source/dataset hash entries referring to the private research layout were removed from the public manifest; its original hash is retained in `paper_provenance.json`. The application manifest retains specification names and fit counts. Original observation data, correspondence, private methods drafts, referee exchanges, and development archives are not included. Literature PDFs are cited rather than redistributed.

## Diagnostic panels

- `results/manuscript_v09/`: candidate/active-constraint audits and interior-only coverage. Interpret coverage conditional on eligible successful interior fits, with the denominators shown in the paper.
- `results/manuscript_v010/`: normal-reference boundary/interior coverage records and summaries, plus aggregate clinical derivative/SE diagnostics. Ambient, fixed-face and rank-tangent calculations remain diagnostics, not a boundary inference theorem.
- `results/manuscript_v011/`: Student-reference comparison using the same fitted estimates and SEs. It changes the multiplier only and retains the normal results for comparison.

Run `python tools/recompute_multiplier.py` to regenerate the multiplier comparison directly from the included simulated interval records. No new histories or fits are needed. Monte Carlo counts, failed-fit denominators, and the distinctions between pooled/interior/boundary strata must be retained in any further analysis.

## Building and checking

```bash
python -m unittest discover -s tests -v
python examples/quickstart.py
python examples/time_dependent.py
python examples/alternating.py
python examples/from_csv.py
python tools/recompute_multiplier.py
python -m pip install '.[docs]'
python tools/build_docs.py
python tools/build_paper.py
```

`tools/build_paper.py` requires Tectonic on `PATH`, compiles from the `paper/` directory, and atomically replaces `paper/paper.pdf` after a successful build. Tectonic resolves the document's LaTeX packages and performs the reruns needed for references. The source deliberately uses no private files: figures and linked supplements resolve inside this repository.

`tools/build_docs.py` generates portable HTML documentation and copies the paper source, PDF, and supplements into `site/`. It uses MathJax 3.2.2 from a CDN by default; `--mathjax-url` selects a locally hosted installation. No website deployment occurs. `paper/web-version.md` is a retained accessible HTML-oriented mirror of manuscript 0.12; it is not the canonical source.

The recorded test environment is in `docs/validation_environment.json`; dependency lower bounds are compatibility targets, not claims that every version has been exercised locally. CI covers additional Python versions when run by the repository host.
