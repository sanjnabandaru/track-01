# [Team Name]: Not Just Born Too Soon

**How much of the early-onset preeclampsia placenta signature is delivery timing, and which cells carry the rest?**

OpenBio 2026 Datathon, Track 1: Maternal-Fetal Cell Communication in Preeclampsia

![Headline: cross-validated low-oxygen gene score for every placenta against delivery week](figures/headline.png)

## Team

| Name | Role |
|---|---|
| [Name] | [Role, e.g. analysis lead] |
| [Name] | [Role, e.g. biology and hypothesis] |
| [Name] | [Role, e.g. article and video] |

## Research Question

Starter questions 6 and 9. In villous placenta from early-onset preeclampsia, (Q6) how much of the gene-activity difference from control placentas is explained by gestational age at delivery, and (Q9) are the genes that remain normally expressed across the placenta or concentrated in specific cell types?

## Data Used

- [x] GSE203507 bulk RNA-seq, villous tissue, 30 pregnancies (required)
- [x] Zeisel/Admati gene browser: LEP, FLT1, HTRA4, PAPPA2, ENG, HSD17B1, HLA-G (required, seven genes)
- [x] Admati/Skarbianskis single-cell pseudobulk files (optional extension; used as a healthy reference and for within-cell-type comparisons)
- [x] Other public data: MSigDB Hallmark gene sets v7.0, via the [GSEApy test-data mirror](https://github.com/zqfang/GSEApy/tree/master/tests/extdata) (checksum-verified in the notebook)

## How to Run

The data package from the datathon is included in `data/`, so the repository runs on its own.

```bash
git clone --branch integrate-track01-analysis https://github.com/sanjnabandaru/track-01.git
cd track-01
python3.11 -m venv .venv && source .venv/bin/activate     # CPython 3.11.8 used for the submitted run
pip install -r requirements.txt
cd notebooks
jupyter nbconvert --to notebook --execute track01_workflow.ipynb --inplace
jupyter nbconvert --to notebook --execute track01_analysis.ipynb --inplace  # optional advanced companion
```

Once these changes are published, open `notebooks/track01_workflow.ipynb` from this branch in Google Colab and run all cells. The first cell clones the published branch and installs `requirements.txt`; unpublished local changes are not available in Colab.

Reproducibility details:

- Core notebook: `notebooks/track01_workflow.ipynb` follows the official workflow and compares 15 PE placentas with 5 spontaneous preterm controls. Its outputs go to `notebooks/results/workflow/`.
- Advanced companion: `notebooks/track01_analysis.ipynb` adds adjusted models and cell-type analyses; `notebooks/track01_analysis.py` contains the same code. Its 15-vs-10 preterm model also includes 5 FGR-only controls and is a distinct, secondary comparison.
- Browser notes read by both notebooks: `notebooks/inputs/track01_browser_lookups.csv`. They are qualitative annotations, not donor-level validation; the team should verify each entry in the browser before submission.
- Environment: every package pinned in `requirements.txt`. Both notebooks were rerun from scratch in a fresh CPython 3.11.8 virtual environment. The main notebook takes about 15 minutes on 8 cores; set `RUN_SLOW = False` to skip the leave-one-out and permutation loop.
- Random seed: 2026 (permutations and null gene sets). Small numerical differences can remain across environments.
- Internet is needed once for Hallmark v7.0 (about 50 KB, checksum-verified, not redistributed here). A local `data/track01_gene_sets.gmt` must match that exact checksum; the official `data/get_gene_sets.py` downloads a different release.

## Files

| Path | Contents |
|---|---|
| `notebooks/track01_analysis.ipynb` | Secondary advanced analysis with outputs |
| `notebooks/track01_analysis.py` | Same code as a script |
| `notebooks/track01_workflow.ipynb` | Rubric-facing core workflow |
| `notebooks/inputs/track01_browser_lookups.csv` | Zeisel/Admati browser notes (seven genes) |
| `notebooks/results/` | All result tables (`workflow/` holds the walkthrough's outputs) |
| `figures/` | Graphs from the main analysis, including `headline.png` |
| `docs/article.md` | Science communication article ([HUDP guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/main/docs/article-guide.md)) |
| `docs/analysis-plan.md` | Comparisons, models, thresholds and confounders |
| `docs/hypothesis-card.md` | Hypothesis and how to test it |
| `docs/presentation-script.md` | Five-minute talk script |
| `docs/rubric-checklist.md` | Where each rubric requirement is met |
| `requirements.txt` | Pinned Python packages |
| `data/`, `TRACK01_GUIDE.md`, `1-setup.md` to `4-submit.md`, `reference/`, `templates/` | The organizers' data package and track materials, unchanged except setup notes |

## Key Findings

- **Core Q6:** the 15-vs-5 PE-versus-spontaneous-preterm two-group fit finds **472 genes** at the stated threshold. A separate three-group fit finds 532 for the same contrast; the fits estimate dispersion using different samples, so their counts should not be mixed. Swapping in term controls also gives 532 genes but only 240 shared; fold changes still correlate at 0.76, with 95% in the same direction.
- **Secondary Q6:** adding timing covariates changes the estimated PE effects by a median **8%** (term indicator) or **2%** (linear weeks) among 266 selected genes. These model-dependent changes are not a causal estimate of the fraction due to timing. The adjusted 15-vs-10 preterm model finds **190 genes**; it is not the core comparison.
- **Exploratory Q9:** among available genes, 70 of 116 PE-up genes peak in a trophoblast class in healthy reference donors. Donor-level pseudobulk analyses suggest a syncytiotrophoblast signal, but small groups and cell mixture limit localization.
- **Hypothesis:** hypoxia-response genes are enriched and the cross-validated score is about 2.5x higher in PE in this dataset. A non-significant delivery-week slope does not establish that timing has no effect; the proposed mechanism needs independent testing.

## Limitations

- Small groups: 15 PE, 10 preterm non-PE and 5 term placentas; single-cell controls 6 (late) and 3 (early).
- Timing modelled as a line and as a term switch; curved effects within 26 to 34 weeks cannot be ruled out.
- Labour, delivery mode and medications are not recorded in the bulk data. Single-cell checks and an FGR-only control model address this only partly.
- Bulk tissue mixes cell types; cell-type assignment relies on another study.
- Tissue taken at delivery shows association, not cause.

## Presentation

Video link: [add your unlisted YouTube link]

## Data Sources

1. Ackerman WE, Buhimschi IA, Buhimschi CS, et al. Molecular subphenotyping using placental transcriptomics in early-onset preeclampsia and fetal growth restriction. NCBI GEO [GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507).
2. Admati I, Skarbianskis N, et al. Placenta PE single-cell RNA-seq dataset. Figshare, [doi:10.6084/m9.figshare.23264102.v1](https://doi.org/10.6084/m9.figshare.23264102.v1) (CC BY 4.0). Paper: Admati I, Skarbianskis N, Hochgerner H, et al. Two distinct molecular faces of preeclampsia revealed by single-cell transcriptomics. *Med* 4(10):687-709 (2023), [doi:10.1016/j.medj.2023.07.005](https://doi.org/10.1016/j.medj.2023.07.005).
3. Zeisel lab placenta PE gene browser. https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html
4. Liberzon A, et al. The Molecular Signatures Database Hallmark gene set collection. *Cell Systems* 1(6):417-425 (2015). [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
5. Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39(9) (2023). [doi:10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547)
