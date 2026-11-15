# Analysis Plan

## Question

How does the choice of term versus spontaneous-preterm controls affect observed gene-expression differences in early-onset preeclampsia? These analyses are exploratory; this document does not claim preregistration.

## Inputs and unit of analysis

- GSE203507 villous counts and metadata from the released Track 1 data package: 30 samples, one per pregnancy.
- Primary comparison: 15 preeclampsia samples versus 5 spontaneous-preterm controls. Five term controls and five FGR-only samples are described separately; the FGR-only group is excluded from the primary comparison.
- Delivery weeks overlap between cases and spontaneous-preterm controls, but samples are not individually matched. Spontaneous-preterm and term delivery circumstances also differ.
- Zeisel/Admati browser: qualitative cell-type context for five genes; independently verify the notebook's manual annotations.

## Methods

1. Check sample IDs, integer nonnegative counts, duplicate IDs and metadata. Retain the released data and report limitations.
2. Retain genes with at least 10 reads in at least 5 samples; examine PCA and delivery weeks.
3. Fit PyDESeq2 with `~q6_group`, contrasting PE with spontaneous-preterm controls. Report genes with Benjamini-Hochberg adjusted p < 0.05 and absolute log2 fold change >= 1.
4. Explore term and preterm control contrasts using a separate three-group model. Its results are not identical to a two-group fit.
5. Explore 20 seeded five-versus-five resamples, threshold sensitivity and covariates for FGR, fetal sex and delivery week, plus exclusion of S01/S12. These are exploratory checks with limited power, not independent validation.
6. Describe browser observations as qualitative context only. Bulk tissue cannot establish within-cell-type changes.

## Limitations

Small groups, overlapping rather than matched weeks, delivery circumstances, possible cell-mixture differences, and tissue collected only at delivery prevent attribution of a causal timing fraction. Control-list overlap and differences in significance are descriptive, not a decomposition into proven disease and timing effects.
