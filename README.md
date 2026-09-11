# A Metal-Poor Tidal Disruption Engine: The Warm Molecular Equilibrium and Its Application to Little Red Dots

Manuscript source, bibliography, figures, and the equilibrium-solver scripts behind
Guillochon & Loeb. Every power law, table entry, and figure in `paper/lrd_engine.tex` is
produced by the commands below. The solver builds on the one released with Guillochon & Loeb
(2026), arXiv:2608.28947 (https://github.com/guillochon/tde-engine).

## Layout

| Path | Contents |
| --- | --- |
| `paper/lrd_engine.pdf` | The compiled paper, built from the source in this repo at the commit shown in the commit message |
| `paper/lrd_engine.tex`, `paper/lrd.bib` | Manuscript (AASTeX 7.0.1 two-column, `pdftex` class option; `aasjournal.bst` included) and ADS-generated bibliography |
| `paper/twindow.pdf`, `paper/fplane.pdf` | Figures 1 and 2 |
| `scripts/losscone.py` | GL26 bridged Cohn-Kulsrud loss-cone flux and gravitationally focused collision rate (unmodified) |
| `scripts/lrd_engine.py` | Reproduces the GL26 CO-cooled belt exactly, then the fixed-200 K and turbulent-thermal-balance H2 variants that fail (Section 4.2) |
| `scripts/lrd_final.py` | Production solve: power laws in (m6, T3, x_H2, zeta), compression floor, temperature window, f_* scan; writes `scripts/lrd_results.json` and `scripts/lrd_scan.json` |
| `scripts/lrd_figs.py` | Draws the two figures from the JSON outputs |
| `scripts/check_numbers.py` | Evaluates the JSON power laws at the fiducial point and prints every number the text quotes |
| `scripts/lrd_estimates.py` | Order-of-magnitude estimates in the text that are not solver power laws: X-ray heating of the belt, belt dust emission, flare variability, [O I] cooling, H2 formation times, zeta propagation |
| `scripts/build_bib.py`, `scripts/ads_validate.py`, `scripts/verify_citations.py` | Bibliography from ADS bibcodes (`scripts/bibcodes.json`) and the abstract dump used to check relevance |
| `scripts/check_abstract.py`, `scripts/check_floats.py` | arXiv abstract length / macro check; figure and table placement check |
| `scripts/arxiv_zip.py` | Bundles the files arXiv needs |

## Environment

* Python 3.11+ with the packages in `requirements.txt` (`numpy`, `scipy`, `sympy`, `matplotlib`, `requests`).
* An ADS API token in `.env` as `ADS_TOKEN=...` (bibliography only).
* A TeX distribution with `aastex` (7.0.1) and `latexmk`; `aasjournal.bst` is included.

## Reproducing the paper

```bash
# 1. Equilibrium solve (a few minutes; sympy + brentq over the loss-cone code)
python scripts/lrd_final.py

# 2. Figures
python scripts/lrd_figs.py

# 3. The numbers quoted in the text, for checking against the manuscript
python scripts/check_numbers.py
python scripts/lrd_estimates.py

# 4. Bibliography from ADS (needs .env) and manuscript
export $(cat .env | tr -d '\r')
python scripts/build_bib.py
python scripts/check_abstract.py
cd paper && latexmk -pdf lrd_engine.tex
```

The JSON outputs from the paper's run are included so steps 2 and 3 work without step 1.
`scripts/lrd_engine.py` is not part of the pipeline; it documents the two failed closures
for the temperature discussed in Section 4.2 and checks that the CO branch reproduces GL26.
