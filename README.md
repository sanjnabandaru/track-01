# Not Just Born Too Soon

**Track 1: How does the choice of control group affect the placental gene-expression differences observed in early-onset preeclampsia?**

## Team

Team: **Reads between the Lines**

Members and roles: **[add every member's name and role]**

Presentation: **[add unlisted YouTube link; every team member must appear or speak]**

## Question and data

Preeclampsia often leads to early delivery, so differences from term placentas could reflect delivery timing or other differences between pregnancies. We compare 15 placentas from pregnancies with early-onset preeclampsia with 5 spontaneous-preterm controls, then explore term controls and sensitivity models. These observational groups have overlapping delivery weeks but are not individually matched.

The released [GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507) villous-tissue data contain 30 samples, one per pregnancy: 15 with preeclampsia, 5 with fetal growth restriction only, 5 spontaneous-preterm controls, and 5 term controls. The notebook records qualitative notes for five genes from the [Zeisel/Admati placental gene browser](https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html); the team must verify those notes before submission.

## Run the analysis

The data needed by the notebook are included in `data/`. The notebook was run with CPython 3.11.8 and the versions in `requirements.txt`.

```bash
git clone https://github.com/sanjnabandaru/track-01.git
cd track-01
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cd notebooks
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=3600 track01_workflow.ipynb
```

Alternatively open `notebooks/track01_workflow.ipynb` in Jupyter or Google Colab and restart and run all cells. Colab clones the repository's `main` branch and installs its pinned requirements. Outputs are saved under `notebooks/results/`; the article's exported images are in `figures/`. Random draws use seed 2026. Runtime and numerical results near the significance cutoff can vary by machine.

## Submission files

| Path | Contents |
|---|---|
| `notebooks/track01_workflow.ipynb` | Primary analysis and outputs |
| `notebooks/results/` | Generated tables and figures |
| `figures/headline.png` | Headline visual for the article and video |
| `docs/article.md` | Science communication article |
| `docs/analysis-plan.md` | Analysis plan |
| `docs/hypothesis-card.md` | Testable hypothesis and alternatives |
| `docs/presentation-script.md` | Draft five-minute presentation |
| `docs/rubric-checklist.md` | Remaining submission checks |
| `data/` and `requirements.txt` | Released input files and pinned software |

## Interpretation and limitations

The main DESeq2 model compares preeclampsia with spontaneous-preterm controls; 472 genes pass at Benjamini-Hochberg adjusted p < 0.05 and absolute log2 fold change >= 1. A separate three-group fit gives 532 for that contrast; it is a different model, not a rerun of the primary fit. Additional models and random draws are exploratory. Control groups are small, spontaneous preterm delivery is not a healthy pregnancy control, and delivery circumstances and cell composition may differ. The analysis cannot determine a causal fraction attributable to timing or assign bulk-expression changes to a particular cell type.

## Before submitting

- Add the members, roles, and unlisted video link above and the article byline.
- Independently verify the browser observations and review the scientific claims.
- Use the same `figures/headline.png` in the article and the recording; every team member must appear or speak.
- Submit on [Devpost](https://openbio-datathon.devpost.com/) or by the [official pull-request route](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/main/submissions/README.md).

## Sources

- Ackerman et al., placental transcriptomics in early-onset preeclampsia and fetal growth restriction, [GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507).
- [OpenBio Datathon 2026 Track 1 data package](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/main/track-01/data).
- [Zeisel/Admati placental gene browser](https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html).
- Muzellec et al., [PyDESeq2](https://doi.org/10.1093/bioinformatics/btad547).
