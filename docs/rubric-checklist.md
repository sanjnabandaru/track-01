# Rubric Checklist (Track 1)

Where each judging criterion and submission requirement is met.

## Judging Criteria

| Criterion | Weight | Where we meet it |
|---|---:|---|
| **Biological question and rationale** | 20% | Two focused starter questions (Q6 timing, Q9 cell types), each tied to preeclampsia: README "Research Question"; article Introduction; analysis plan "Research Question"; notebook Section 1 |
| **Data rigor and reproducibility** | 25% | Integrity checks and audit flags reported, not corrected; donor as unit; methods and exploratory analyses documented without claiming preregistration; seed 2026; pinned Python dependencies; checksum-verified v7.0 gene sets; single notebook executed from scratch |
| **Analytical quality** | 25% | Delivery timing handled four ways (PCA audit, model comparison with two timing codings, variance partition, preterm comparison with overlapping delivery weeks) plus a healthy-only timing signature with a noise and power baseline and a worst-case term-only decomposition; FGR and fetal sex in every adjusted model; threshold grid; sensitivity models incl. labour proxy; leave-one-out; label permutations; matched null for specificity; within-cell-type null; single-cell checks for delivery mode and magnesium sulfate; cross-validated headline score (Sections 3 to 11) |
| **Hypothesis quality** | 20% | `docs/hypothesis-card.md`: specific result, mechanism, five alternative explanations, four falsifiable follow-ups, confidence level split by claim, "suggests" language |
| **Communication and visual clarity** | 10% | 12 figures with plain-language titles, axis labels, captions and alt text; article written for general readers with every term explained; colour-blind-safe palette plus marker shapes; five-minute script with timings |

## Submission Requirements

| Requirement | Location | Status |
|---|---|---|
| Reproducible code that runs top to bottom, with a README | `notebooks/track01_analysis.ipynb` (single notebook), `README.md`, `requirements.txt` | Done, rerun from scratch in CPython 3.11.8 |
| Recorded presentation, 5 minutes or less, unlisted YouTube | Script in `docs/presentation-script.md`; link goes in README | Team to record |
| Article following the HUDP guide, at least four graphs | `docs/article.md` (10 graphs) | Draft; team review required |
| Headline visual used in video and article | `figures/headline.png` | Done |
| Completed analysis plan | `docs/analysis-plan.md` | Draft; team review required |
| One completed hypothesis card | `docs/hypothesis-card.md` | Draft; team review required |
| Gene browser, at least three genes | `notebooks/inputs/track01_browser_lookups.csv` (seven annotations) | Verify observations and document access before submission |
| Every dataset cited | README "Data Sources"; article "Data and Code" | Done |

## Before Submitting

- [ ] Team names and roles in README and article byline
- [ ] Independently verify the browser observations and review all scientific claims
- [ ] Video recorded, every member speaks, uploaded unlisted, link in README and opens in a private window
- [ ] Repository link added to article "Data and Code"
- [ ] Submitted on Devpost or by pull request before Sun 15 Nov, 6:00 PM EST
