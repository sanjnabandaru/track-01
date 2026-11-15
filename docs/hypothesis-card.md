# Hypothesis Card

## Title

Early-onset preeclampsia may involve a trophoblast-associated hypoxia-response program, potentially driven by HIF in syncytiotrophoblast.

## Clinical/Biological Problem

Early-onset preeclampsia forces delivery before 34 weeks, so placentas from affected pregnancies are also younger than most controls. Studies that compare them with term placentas cannot tell disease genes from pregnancy-stage genes (starter question 6). The genes that remain also need a cellular home (starter question 9). The placenta releases factors such as soluble FLT1 (sFlt-1) and soluble endoglin (ENG) that damage the mother's blood vessels ([Maynard et al. 2003](https://doi.org/10.1172/JCI17189); [Venkatesha et al. 2006](https://doi.org/10.1038/nm1429)), and the sFlt-1/PlGF ratio is used clinically ([Zeisler et al. 2016](https://doi.org/10.1056/NEJMoa1414838)). Knowing whether these signals are disease or timing, and which cells make them, decides where to look for mechanisms and biomarkers.

## Dataset(s) and Relevant Cells

- GSE203507 villous placenta, bulk RNA-seq, 30 pregnancies: 15 early-onset PE (7 with growth restriction), 5 growth restriction without PE, 5 spontaneous preterm, 5 term.
- Admati/Skarbianskis single-cell pseudobulk, 26 donors: 9 control donors as a healthy reference; late-onset arm (7 PE vs 6 controls, 37 to 40 weeks); early-onset arm (10 vs 3, exploratory).
- Zeisel/Admati gene browser: LEP, FLT1, HTRA4, PAPPA2, ENG, HSD17B1, HLA-G.
- Cells: syncytiotrophoblast (STB) and extravillous trophoblast (EVT).

## Observed Result

**Q6, timing.**
- PE placentas were delivered 3.0 weeks earlier on average, yet the main axis of variation (PC1, 22%) tracks PE (p = 4e-5) and not delivery week (p = 0.14).
- Across 266 PE genes, adding delivery timing to the model removes a median **8%** of each gene's difference (term vs preterm indicator) or **2%** (weeks as a line). Only 1 gene loses more than half. All keep their direction.
- In those genes, PE explains a median 26% of variation and delivery week 1% (partial R squared).
- The core 15-vs-5 PE-vs-spontaneous-preterm two-group fit finds **472 genes**; a distinct three-group fit finds 532 for the same contrast. The secondary adjusted model compares 15 PE with 10 non-PE preterm placentas (including 5 FGR-only) and finds **190 genes** (156 higher, 34 lower). A term-vs-spontaneous-preterm control contrast yields 90 genes, of which 4 overlap; it is not a pure timing effect.
- An exploratory decomposition of the PE-vs-term contrast gives a median timing-related component of 25%; this estimate depends on model assumptions and is not causal.

**Q9, cell types.**
- In healthy placenta, the genes that rise in PE are **concentrated, not broad**: median specificity (tau) 0.69 vs 0.31 for expression-matched random genes. 70 of 116 peak in a trophoblast cell type, vs 30% of all genes. Timing genes are also specific, but not trophoblast-enriched.
- In the single-cell dataset, the bulk "higher" genes show the strongest median increase in **syncytiotrophoblast**: 81% are higher in late-onset PE donors with overlapping delivery weeks (median +0.47 log2), and 96% in early-onset donors. The small-donor subgroup checks are exploratory and do not rule out confounding.
- Bulk STB marker scores do not show a significant PE effect. This does not exclude differences in cell mixture.

**Key genes (95% CI).** FLT1 5.6x (3.0 to 10.3), ENG 3.9x (2.7 to 5.4), LEP 26x (9 to 79), HTRA4 5.3x (3.0 to 9.5), FSTL3 6.1x (2.9 to 13.0), in the secondary preterm model. Adding timing changes the estimated effects by 6% to 26%; this is not a causal fraction.

**Program.** The higher genes are enriched for Hallmark hypoxia (19 genes, adjusted p = 6e-13) and glycolysis. A cross-validated hypoxia score, with genes chosen without the placenta being scored, is about **2.5 times higher in PE** (95% CI 1.7 to 3.8). Its delivery-week slope was not significant; that does not prove a lack of timing effect. The score needs independent validation.

## Proposed Mechanism

We hypothesize that poor remodelling of maternal spiral arteries reduces placental perfusion and stabilizes HIF in syncytiotrophoblast, inducing stress-response genes and secreted factors that contribute to endothelial dysfunction. The observational results are consistent with this mechanism but cannot establish it, determine direction of causation, or exclude residual timing and cell-mixture effects.

## Why the Result Is Biologically Plausible

- Low oxygen raises sFlt-1 in human villous explants, and HIF-1α knockdown lowers it ([Nevo et al. 2006](https://doi.org/10.1152/ajpregu.00794.2005)).
- Hypoxia raises leptin in trophoblast cells, and placental leptin is raised in severe PE ([Mise et al. 1998](https://doi.org/10.1210/jc.83.9.3225)).
- The source single-cell study reported cell-autonomous dysregulation of FLT1 in the syncytium in early-onset PE ([Admati et al. 2023, Med](https://doi.org/10.1016/j.medj.2023.07.005)).
- STB lines the villous surface bathed in maternal blood, so its secreted products reach the mother directly.

## Alternative Explanations/Confounders

1. **Residual timing effects.** Our timing models assume either a straight line in weeks or a single term vs preterm step. A curved timing effect between 26 and 34 weeks could hide in the preterm range. The timing-matched comparison and the flat headline trend argue against a large one.
2. **Labour and delivery mode.** Not recorded in bulk. Spontaneous preterm controls are likely laboured; PE deliveries are mostly pre-labour caesareans. Comparing PE with growth-restriction controls (also mostly caesarean) keeps every gene's direction (100%), with effects about 70% as large (median ratio 0.71), but loses power (5 controls, 11 genes pass). The single-cell STB result holds in caesarean-only donors.
3. **Magnesium sulfate and antihypertensives.** Not recorded in bulk. In single-cell STB, the result holds without magnesium sulfate (4 PE vs 6 controls).
4. **Cell mix.** STB and EVT marker scores do not rise. Many genes peak per cell in EVT, and the healthy EVT profile rests on 3 control donors, so an EVT contribution cannot be excluded.
5. **Reverse direction.** Maternal hypertension could reduce placental perfusion and cause the program, rather than the reverse. Tissue taken at delivery cannot separate these.

## Proposed Validation Experiment or Follow-up Analysis

1. **Computational:** test the same cross-validated hypoxia score in a larger public cohort with delivery mode and labour recorded (for example GSE75010, 157 placentas). Prediction: the PE gap persists after adjusting for gestational age, delivery mode and labour. If it disappears, the hypothesis fails.
2. **Timing curve:** in a cohort with many preterm controls, fit a flexible (spline) gestational-age curve to test for curved timing effects within 26 to 34 weeks.
3. **Spatial:** RNA in situ hybridization for FLT1, LEP and NDRG1 on PE and timing-matched preterm sections, co-stained for STB (CGB) and EVT (HLA-G). Prediction: signal in STB, stronger near hypoxic or infarcted villi.
4. **Mechanistic:** primary trophoblast differentiated to syncytium at 1%, 8% and 20% oxygen, with and without HIF-1α/HIF-2α knockdown. Prediction: low oxygen induces the program and knockdown blocks FLT1, LEP and ENG. If knockdown does not block it, the HIF link is wrong.

## Potential Intervention Relevance

The program includes proteins already measurable in blood (sFlt-1, soluble endoglin, leptin, activin A from INHBA). A panel tied to a trophoblast hypoxia program is a realistic diagnostic follow-up. HIF signalling is druggable, but our data come from delivery and cannot support treatment claims.

## Confidence Level and Limitations

**Medium.**

- **Q6:** the PE comparison shows a signal with preterm controls, but the timing contribution remains uncertain and model-dependent; the analyses are not independent.
- **Moderate confidence (Q9):** the genes are trophoblast-concentrated in healthy placenta and rise most inside STB, including in timing-matched donors. Single-cell control groups are small (6 late, 3 early), and within-STB tests mostly do not reach adjusted significance gene by gene.
- **Lower confidence:** that HIF drives the program and that it causes rather than follows the maternal disease.
- **What expression data cannot establish:** protein secretion, causation, or the early events that start preeclampsia. Bulk samples lack labour, delivery mode and medication records. Groups are small (15 PE, 10 preterm non-PE, 5 term).
