# Reads between the Lines: Presentation Script (five minutes maximum)

Assign parts so **every team member appears or speaks**. Rehearse with the actual team and check the recording length. Use `figures/headline.png` both here and in the article.

## 0:00–0:45 — Question

Preeclampsia often leads to early delivery. When we compare placentas, differences could reflect the disease, delivery timing, or other features of the pregnancies. We asked how the choice of preterm versus term controls changes the differences we observe.

## 0:45–1:45 — Data and methods

We used the released GSE203507 villous-tissue data: 15 early-onset preeclampsia placentas and 5 spontaneous-preterm controls in the primary comparison, plus other control groups for exploratory analyses. Each pregnancy counts once. The preterm groups have overlapping delivery weeks, but are not individually matched. We used PyDESeq2 and required a twofold difference and adjusted p below 0.05.

## 1:45–3:30 — Results

Show the headline image, delivery-week graph, PCA and volcano plot. Explain what each axis and point means. The primary comparison finds gene-expression differences against preterm controls. Other comparisons and seeded five-versus-five draws show how control choice affects results, but do not isolate timing from other delivery circumstances. Explain the numerical results directly from the final notebook.

## 3:30–4:45 — Hypothesis and limits

Our hypothesis is that some differences, including higher FLT1 RNA, will replicate in a larger study with well-characterized preterm controls. This dataset has few controls, bulk tissue mixes cell types, and samples were collected at delivery. We cannot infer cause or a causal fraction due to timing. A follow-up should record gestational age, labour and delivery mode and measure RNA and protein separately.

## Final slide

List every team member, data sources, repository link and acknowledgements. Verify browser observations before presenting them.
