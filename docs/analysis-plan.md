# Analysis Plan

Team: [Team Name]. Track 1. This document records the analyses and their outcomes; it is not a preregistration. The core workflow and the exploratory companion use different comparisons, as described below.

## Research Question

We answer two linked starter questions.

- **Q6:** How much of the difference between early-onset preeclampsia (PE) and control placentas could be explained by when delivery happened?
- **Q9:** Do the genes that differ in preeclampsia normally work across the placenta, or are they concentrated in specific cell types?

Q6 comes first because PE forces early delivery: any PE vs control gene list mixes disease with pregnancy stage. Only the part of the signal that survives a timing check is worth localizing to cells in Q9.

## Datasets

| Dataset ID | Role | Why This Dataset |
|---|---|---|
| GSE203507 villous bulk RNA-seq (`track01_counts_vt.csv`, `track01_metadata_vt.csv`) | Q6 and the gene list for Q9 | 30 placentas, one per pregnancy. The groups overlap in delivery week; the 5 term and 5 spontaneous preterm controls provide a timing-related contrast, not an isolated timing effect |
| `qc/track01_htseq_summary.csv` | Library quality | Read assignment rate for a sensitivity analysis |
| `interpretation/track01_marker_matrix.csv` | Q9 cell-mix check | Marker genes for nine broad cell classes |
| Admati/Skarbianskis pseudobulk (`extension/pseudobulk/`) | Q9 | Control donors give a healthy reference for cell-type specificity. The late-onset arm (7 PE vs 6 controls, 37 to 40 weeks) is matched on delivery timing by design. Donor table records delivery mode and magnesium sulfate |
| Zeisel/Admati gene browser | Q9, qualitative | Required; seven genes looked up |
| MSigDB Hallmark gene sets (v7.0, checksum-verified) | Pathway context and headline score | Predefined sets we did not choose |

We did not use the decidua extension: it comes from the same 27 pregnancies, so it adds no independent donors.

## Methods

1. **Data loading and checks.** Confirm 30 samples in both tables, integer non-negative counts, no duplicates, one sample per donor. Summarize clinical fields. Keep `NA` as `NA` and list `qc_flags` without changing values. Compute read assignment rates.
   - **Outcome:** all checks pass. S12 and S22 have fewer than 25% of reads assigned to genes. PE placentas were delivered 3.0 weeks earlier on average (30.2 vs 33.3). No PE delivery after 33.7 weeks; no term birth before 39.0 weeks. Fetal sex is unbalanced: 10 of 15 PE vs 5 of 15 non-PE female.

2. **Exploration.** 15,835 genes with at least 10 reads in at least 5 samples, log2 CPM, PCA. Test PC1 to PC5 against every clinical field instead of judging colours by eye.
   - **Outcome:** PC1 (22%) tracks PE (p = 4e-5), not delivery week (p = 0.14). Delivery week lines up with PC2 (p = 0.01), with growth restriction and sequencing depth.

3. **Q6, model comparison on all 30 placentas.** PyDESeq2 0.5.4, four designs, PE vs no PE: naive `~ pe_status`; no timing `~ fetal_sex + fgr_status + pe_status`; plus weeks at delivery (linear); plus a term-birth indicator. For each gene passing the no-timing model, compute 1 minus (log2 fold change with timing / log2 fold change without timing). This selected-gene, model-dependent contrast is not a causal fraction attributable to timing.
   - **Outcome:** median change is 8% with the term indicator and 2% with linear weeks. 1 of 266 genes loses more than half its estimated difference. Every gene keeps its direction.

4. **Q6, variance partition.** Per-gene linear model of log2 CPM on PE, FGR, fetal sex and weeks (all 30). Partial R squared for each factor.
   - **Outcome:** in the 266 PE genes, PE explains a median 26% of variation and weeks 1%.

5. **Q6, core and secondary comparisons.** The rubric-facing workflow compares 15 PE placentas with 5 spontaneous preterm controls; its two-group fit finds 472 passing genes, while its separate three-group fit finds 532 for the same contrast. The exploratory companion compares PE vs 10 preterm non-PE placentas (5 spontaneous preterm and 5 FGR-only), with `~ fetal_sex + gestational_age + fgr_status + pe_status`. The groups overlap in delivery weeks but are not individually matched. A gene passes at BH adjusted p < 0.05 **and** |log2 fold change| >= 1. The term-vs-spontaneous-preterm control contrast (5 vs 5) also differs in delivery circumstances.
   - **Outcome:** 190 genes (156 higher, 34 lower). Timing signature: 90 genes, of which 4 overlap the PE list. Fold changes agree with the all-30 no-timing model (r = 0.94).
   - **Added after comparing with a teammate's parallel analysis:** a noise and power baseline (8 random 5-vs-5 splits of the controls, and 8 draws of 5 PE vs the 5 spontaneous preterm births) and a worst-case term-only decomposition (three-group model, PE vs term split into PE vs preterm plus preterm vs term). **Outcome:** noise 0 to 18 genes, disease 185 to 420, timing 90; term-only median timing part 25%.

6. **Robustness.** Threshold grid; sensitivity models (a) drop flagged or low-QC samples, (b) spontaneous preterm controls only, (c) add maternal age and depth, (d) growth-restriction controls only as a proxy for matching labour and delivery mode; leave-one-out over the 10 preterm controls; 10 label permutations within FGR strata (seed 2026).

7. **Q9, cell mix in bulk.** Marker scores per class, PE effect with the same covariates (linear regression, 95% CI).

8. **Q9, specificity in healthy placenta.** Pseudobulk control donors only, columns with at least 10 cells. Tau specificity index per gene (0 = equal in all cell types, 1 = one cell type). Compare PE genes with 1,000 expression-matched random sets, and with the delivery-timing genes. Count which class each gene peaks in.
   - **Outcome:** PE "higher" genes median tau 0.69 vs 0.31 for matched random genes; 70 of 116 peak in trophoblast vs 30% of all genes. Timing genes are also specific but not trophoblast-enriched.

9. **Q9, which cell type changes.** Within each class, PE minus control mean log2 CPM over the bulk "higher" genes, donor as unit, against 500 same-size random sets. Late onset (timing-matched) and early onset (exploratory). In late-onset syncytiotrophoblast, repeat for caesarean-only donors and for donors without magnesium sulfate, and fit DESeq2 with delivery mode as a covariate.

10. **Gene browser.** LEP, FLT1, HTRA4, PAPPA2, ENG (higher), HSD17B1 (lower) and HLA-G (a delivery-timing gene). Chosen for high baseMean, stability in leave-one-out and to represent the story.

11. **Headline.** Cross-validated hypoxia score: for each preterm placenta, Hallmark hypoxia genes are chosen on the other 24, then the held-out placenta is scored. Term births are scored with genes chosen from all 25. Plotted against delivery week.

## Confounders and Controls

| Confounder | Why it matters | How we handle it |
|---|---|---|
| Gestational age at delivery | The question itself (Q6) | Four-model comparison with two timing parametrizations; variance partition; timing-matched preterm-only model; separate timing signature from healthy controls; headline plotted against delivery week |
| Fetal growth restriction | 7 of 15 PE also FGR; FGR lines up with PC2 | Covariate in every adjusted model; FGR-only control model |
| Fetal sex | 67% female in PE vs 33% in non-PE | Covariate |
| Sequencing depth, library quality | Depth lines up with PC2; two low-assignment samples | DESeq2 size factors; sensitivity models (a) and (c) |
| Maternal age | 15 to 39 years | Sensitivity model (c) |
| Labour and delivery mode | Not recorded in bulk; spontaneous preterm births are usually laboured | Bulk: FGR-only control model. Single-cell: caesarean-only subset and delivery-mode covariate |
| Magnesium sulfate | Given to most PE patients; not recorded in bulk | Single-cell: subset without magnesium sulfate |
| Audit flags | S04, S15 reached 140/90 without PE; S12 PE never did | Kept, reported, dropped in sensitivity (a) |
| Cell mix | Bulk tissue mixes cell types | Marker scores; within-cell-type comparison |
| Donor non-independence | Cells and samples from one person are not independent | One villous sample per pregnancy; donor x class pseudobulk, never cells |

## Expected Figures

1. Headline: cross-validated hypoxia score for every placenta against delivery week (`figures/headline.png`)
2. PCA by group and by delivery week (`fig2_pca.png`); PC-by-metadata audit (`fig3_pc_audit.png`)
3. Q6: fold changes with vs without timing and share explained (`fig4_timing_share.png`); variance partition (`fig5_variance_partition.png`)
4. Timing-matched volcano (`fig6_volcano.png`); Hallmark enrichment (`fig7_pathways.png`); robustness (`fig8_robustness.png`)
5. Q9: specificity and peak cell type (`fig9_specificity.png`); within-cell-type change (`fig10_within_celltype.png`); key-gene heat map (`fig11_celltype_heatmap.png`)
6. FLT1, ENG and LEP in every placenta (`fig12_key_genes.png`)
