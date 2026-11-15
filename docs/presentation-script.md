# Presentation Script (5 minutes)

Four speakers (A to D); merge parts for a smaller team. Every member must appear or speak at least once. Rehearse to about 4:45. Show `figures/headline.png` at the start and the end.

## Slide 1, Introduction (0:00 to 0:45). Speaker A

**On screen:** title, headline graph.

> Preeclampsia is dangerously high blood pressure in pregnancy, and it is thought to start in the placenta. In its early form, the only cure is delivery, often of a very premature baby.
>
> That creates a trap for anyone studying these placentas: they are also younger. So we asked two of the starter questions. Question 6: how much of the preeclampsia signal is just delivery timing? Question 9: are the genes that remain spread across the placenta, or concentrated in particular cells?
>
> We used 30 placentas from GSE203507, single-cell data from 26 donors, and the Zeisel gene browser.

## Slide 2, Methodology (0:45 to 1:45). Speaker B

**On screen:** a timeline of delivery weeks by group (or `fig3_pc_audit.png`), then a sketch: "fold change without timing vs with timing".

> Preeclampsia placentas were delivered three weeks earlier on average, and none after 34 weeks.
>
> Our core comparison is 15 preeclampsia placentas against 5 spontaneous preterm controls. Their delivery weeks overlap, but are not individually matched. We also explored PCA and adjusted models, including a broader preterm comparison that adds 5 growth-restriction-only controls.
>
> The advanced models adjust for growth restriction and fetal sex; the core fit is unadjusted. A gene counts if it changes at least 2-fold with adjusted p below 0.05. Each pregnancy is one sample. These analyses are exploratory, not preregistered.

## Slide 3, Results for Question 6 (1:45 to 2:45). Speaker C

**On screen:** `fig2_pca.png`, then `fig4_timing_share.png`, then `fig6_volcano.png`.

> The core two-group fit finds 472 passing genes. A separate three-group fit finds 532 for the same contrast; the two fits use different samples to estimate variation.
>
> In a secondary analysis of 266 selected genes, adding timing changed the estimated effect by a median of 8 percent with a term indicator, or 2 percent with linear weeks. These are model-dependent summaries, not causal estimates.
>
> The broader adjusted preterm model finds 190 genes. A separate term-versus-preterm control contrast finds 90, with 4 shared. Small groups and differences in delivery circumstances mean we cannot isolate a pure timing effect.

## Slide 4, Results for Question 9 (2:45 to 3:45). Speaker D

**On screen:** `fig9_specificity.png`, `fig10_within_celltype.png`, a browser screenshot of LEP labelled "qualitative context", then `fig7_pathways.png`.

> In healthy placenta, these genes are concentrated, not broad: specificity 0.69 versus 0.31 for matched random genes. 70 of 116 peak in trophoblast, against 30 percent of all genes.
>
> Donor-level analyses suggest a signal in syncytiotrophoblast, the cell layer bathed in the mother's blood. The donor groups are small and have overlapping, not individually matched, delivery weeks.
>
> Hypoxia-response genes are enriched. Our cross-validated headline score is about 2.5 times higher in preeclampsia in this dataset, but needs testing in an independent cohort.

## Slide 5, Discussion (3:45 to 5:00). Speaker A, or each member one line

**On screen:** hypothesis in one sentence; three limitations; three next steps.

> Our hypothesis is that, in early preeclampsia, syncytiotrophoblast responds to low oxygen through HIF, raising FLT1, ENG and leptin. We cannot yet establish that mechanism or exclude delivery timing.
>
> Limits: 15 versus 10 placentas; timing modelled as a line and a switch; labour and medications not recorded in the main study; association, not cause.
>
> To test it: rerun the score in a larger cohort that records delivery mode; stain tissue for FLT1 and leptin with trophoblast markers; and grow placental cells at low oxygen with HIF blocked. If blocking HIF does not stop FLT1 and leptin from rising, we are wrong.

## Final Slide

Team names; data sources (GSE203507; Admati, Skarbianskis et al., Figshare 23264102; Zeisel lab gene browser; MSigDB Hallmark); repository link.
