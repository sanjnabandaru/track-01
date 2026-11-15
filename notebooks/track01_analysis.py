# %% [markdown]
# # How Much of the Preeclampsia Placenta Signal Is Delivery Timing, and Which Cells Carry the Rest?
#
# OpenBio 2026 Datathon, Track 1. This notebook runs top to bottom from the released data package (`track-01/data/`).
# The article combines the core workflow with exploratory results produced here.
#
# **Starter question 6.** How much of the difference between preeclampsia (PE) and control placentas could be explained by
# when delivery happened?
#
# **Starter question 9.** Do the genes that differ in preeclampsia normally work across the placenta, or are they
# concentrated in specific cell types?
#
# **Why these two together.** Early-onset PE forces early delivery, so every PE vs control comparison risks confusing
# disease with pregnancy stage (Q6). Only the part of the signal that survives that check is worth localizing to cells (Q9).
#
# | Part | What it answers | Sections |
# |---|---|---|
# | A | Data checks and exploration | 2, 3 |
# | B | Q6: how much is delivery timing? | 4, 5, 6 |
# | C | Robustness of the preterm signal | 7 |
# | D | Q9: broad or cell-type concentrated? | 8, 9, 10 |
# | E | Headline and hypothesis | 11, 12 |
#
# Runtime: about 15 minutes on a laptop. Set `RUN_SLOW = False` to skip the leave-one-out and permutation loop (about 8 minutes).
#
# This exploratory companion complements `track01_workflow.ipynb`, the rubric-facing core analysis. Its adjusted
# preterm model uses additional FGR-only controls and should not be treated as the same comparison as the core model.
# Setup: Python 3.11.8 and `pip install -r ../requirements.txt`.

# %% [markdown]
# ## 0. Setup

# %%
import os, sys, subprocess, hashlib, warnings
from pathlib import Path

if "google.colab" in sys.modules:
    if not Path("/content/track-01").exists():
        subprocess.run(["git", "clone", "-q", "--branch", "integrate-track01-analysis",
                        "https://github.com/sanjnabandaru/track-01.git", "/content/track-01"], check=True)
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "-r", "/content/track-01/requirements.txt"], check=True)
    os.chdir("/content/track-01/notebooks")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
SEED = 2026
rng = np.random.default_rng(SEED)
RUN_SLOW = True

# The notebook expects the released data package at ../data (this repository's layout).
# Also found automatically from track-01/<team>/notebooks or submissions/<team>/notebooks. Override with TRACK01_DATA.
DATA_DIR = Path(os.environ.get("TRACK01_DATA", "../data"))
for candidate in [DATA_DIR, Path("../../data"), Path("../../../track-01/data")]:
    if (candidate / "track01_counts_vt.csv").exists():
        DATA_DIR = candidate
        break
assert (DATA_DIR / "track01_counts_vt.csv").exists(), "Data package not found: set TRACK01_DATA to the track-01/data folder"
RESULTS_DIR = Path("results"); RESULTS_DIR.mkdir(exist_ok=True)
FIG_DIR = Path("../figures"); FIG_DIR.mkdir(exist_ok=True)
INF = DefaultInference(n_cpus=os.cpu_count())

# Three categorical colours (colour-blind-safe as a set), with marker shape as a second cue
C_PE, C_PRE, C_TERM = "#eb6834", "#2a78d6", "#1baf7a"
INK, MUTED, GREY = "#0b0b0b", "#52514e", "#c9c8c3"
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titleweight": "bold",
                     "axes.titlesize": 13, "legend.frameon": False})
GROUP_STYLE = {
    "EOPE":     dict(color=C_PE,   marker="o", label="Preeclampsia only"),
    "EOPE+FGR": dict(color=C_PE,   marker="^", label="Preeclampsia + growth restriction"),
    "FGR":      dict(color=C_PRE,  marker="^", label="Growth restriction, no preeclampsia"),
    "iPTB":     dict(color=C_PRE,  marker="o", label="Spontaneous preterm birth"),
    "TB":       dict(color=C_TERM, marker="s", label="Term birth"),
}
def save(fig, name):
    fig.savefig(FIG_DIR / name, bbox_inches="tight", facecolor="white")
    plt.show()
def nice(cls):
    return cls.replace("_", " ")
print("Data folder:", DATA_DIR)

# %% [markdown]
# ## 1. Questions and Analysis Decisions
#
# Analysis choices and outcomes are summarized in `docs/analysis-plan.md`; this is not a preregistration.
#
# * **Q6 approach.** Estimate the PE effect for every gene twice on all 30 placentas: once *without* delivery timing in
#   the model and once *with* it. The relative change in the estimated log2 fold change is model-dependent, not a causal
#   fraction explained by timing. We model timing two ways, as weeks at delivery (linear) and as a term-birth
#   indicator (which lets the 5 term births differ in any way at all), so the answer does not hang on one assumption.
#   Both models also adjust for growth restriction (FGR) and fetal sex. A third, design-based check compares PE only with
#   preterm placentas delivered in overlapping weeks (26 to 34).
# * **Q9 approach.** Take the genes that survive the timing checks and ask, in *healthy* (control) donors of an independent
#   single-cell dataset, whether each gene is made by many cell types or concentrated in one (the tau specificity index,
#   0 = equal in all cell types, 1 = only one). Compare with expression-matched random genes. Then ask in which cell type
#   the genes change with PE, using late-onset donors with overlapping delivery weeks.
# * **Thresholds.** Benjamini-Hochberg adjusted p < 0.05 and |log2 fold change| >= 1 (at least 2-fold).
# * **Unit of replication.** The donor. One villous sample per pregnancy; in single-cell data each donor x cell class is one
#   pseudobulk column, never individual cells.

# %% [markdown]
# ## 2. Load and Check the Data

# %%
counts = pd.read_csv(DATA_DIR / "track01_counts_vt.csv", index_col=0)
metadata = pd.read_csv(DATA_DIR / "track01_metadata_vt.csv", index_col=0)
htseq = pd.read_csv(DATA_DIR / "qc" / "track01_htseq_summary.csv", index_col=0)
markers = pd.read_csv(DATA_DIR / "interpretation" / "track01_marker_matrix.csv")

checks = {
    "samples in both tables": len(set(counts.columns) & set(metadata.index)),
    "counts nonnegative": bool((counts.values >= 0).all()),
    "counts whole numbers": bool((counts.values % 1 == 0).all()),
    "duplicate genes": int(counts.index.duplicated().sum()),
    "duplicate samples": int(counts.columns.duplicated().sum()),
    "one sample per donor": bool(metadata["donor_id"].is_unique),
}
for k, v in checks.items():
    print(f"{k:28s} {v}")
assert checks["samples in both tables"] == 30 and checks["one sample per donor"]
metadata = metadata.loc[counts.columns]

lib = htseq.loc[counts.columns].copy()
lib["assigned_M"] = lib["assigned_to_genes"] / 1e6
lib["frac_assigned"] = lib["assigned_to_genes"] / lib[["assigned_to_genes", "no_feature", "ambiguous", "alignment_not_unique"]].sum(axis=1)
metadata["library_M"] = lib["assigned_M"]; metadata["frac_assigned"] = lib["frac_assigned"]
low_qc = lib.index[lib["frac_assigned"] < 0.25].tolist()
print("\nReads assigned to genes (millions): median %.1f, range %.1f to %.1f" %
      (lib.assigned_M.median(), lib.assigned_M.min(), lib.assigned_M.max()))
print("Samples with < 25% of reads assigned to genes (kept; dropped in a sensitivity model):", low_qc)

# %%
overview = metadata.groupby("group").agg(
    n=("donor_id", "size"), ga_min=("ga_delivery_weeks", "min"), ga_median=("ga_delivery_weeks", "median"),
    ga_max=("ga_delivery_weeks", "max"), female=("fetal_sex", lambda s: int((s == "female").sum())),
    maternal_age_median=("maternal_age_years", "median"), sbp_max_median=("sbp_max_mmhg", "median"))
display(overview)
display(metadata.loc[metadata.qc_flags != "none", ["group", "qc_flags"]])
flagged_bp = metadata.index[metadata.qc_flags.str.contains("pe_without|no_pe_with")].tolist()
ga_pe = metadata.query("pe_status=='PE'").ga_delivery_weeks.mean(); ga_np = metadata.query("pe_status=='no_PE'").ga_delivery_weeks.mean()
latest_pe = metadata.loc[metadata.pe_status == "PE", "ga_delivery_weeks"].max()
earliest_tb = metadata.loc[metadata.group == "TB", "ga_delivery_weeks"].min()
print(f"Mean delivery week: PE {ga_pe:.1f}, all non-PE {ga_np:.1f} (gap {ga_np-ga_pe:.1f} weeks). "
      f"Latest PE delivery {latest_pe:.1f} weeks; earliest term birth {earliest_tb:.1f} weeks.")
print("Fetal sex: %d of 15 PE pregnancies carried a female fetus vs %d of 15 non-PE." %
      ((metadata.query("pe_status=='PE'").fetal_sex == "female").sum(), (metadata.query("pe_status=='no_PE'").fetal_sex == "female").sum()))

# %% [markdown]
# ## 3. Explore the Samples
#
# Genes with at least 10 reads in at least 5 samples, log2 counts per million (CPM), PCA. Instead of judging colours by
# eye, every principal component is tested against every clinical field (Spearman for numbers, Kruskal-Wallis for groups).

# %%
expressed = (counts >= 10).sum(axis=1) >= 5
log_cpm_all = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
log_cpm = log_cpm_all[expressed]
print(f"{expressed.sum():,} of {len(counts):,} genes pass the expression filter")
pca = PCA(n_components=5, random_state=SEED)
scores = pd.DataFrame(pca.fit_transform(log_cpm.T), index=counts.columns, columns=[f"PC{i}" for i in range(1, 6)])
var = pca.explained_variance_ratio_

def assoc_p(pc, field):
    x = metadata[field]; ok = x.notna()
    if pd.api.types.is_numeric_dtype(x):
        return stats.spearmanr(scores.loc[ok, pc], x[ok]).pvalue
    return stats.kruskal(*[scores.loc[ok & (x == l), pc] for l in x[ok].unique()]).pvalue

fields = {"pe_status": "Preeclampsia", "ga_delivery_weeks": "Weeks at delivery", "fgr_status": "Growth restriction",
          "fetal_sex": "Fetal sex", "maternal_age_years": "Maternal age", "library_M": "Sequencing depth",
          "frac_assigned": "Read assignment rate"}
assoc = pd.DataFrame({pc: [assoc_p(pc, f) for f in fields] for pc in scores.columns}, index=fields.values())
assoc.to_csv(RESULTS_DIR / "pc_metadata_association.csv")
display(assoc.map(lambda p: f"{p:.1e}"))
pe_only_samples = metadata.index[metadata.pe_status == "PE"]
print("Within the 15 PE placentas, PC1 vs delivery week: rho = %.2f, p = %.2f" %
      stats.spearmanr(scores.loc[pe_only_samples, "PC1"], metadata.loc[pe_only_samples, "ga_delivery_weeks"]))

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axes[0]
for g, st in GROUP_STYLE.items():
    idx = metadata.index[metadata.group == g]
    ax.scatter(scores.loc[idx, "PC1"], scores.loc[idx, "PC2"], s=70, edgecolor="white", linewidth=1.2,
               facecolor=st["color"], marker=st["marker"], label=f'{st["label"]} (n={len(idx)})')
ax.set_xlabel(f"PC1 ({var[0]:.0%} of variation)"); ax.set_ylabel(f"PC2 ({var[1]:.0%} of variation)")
ax.set_title("PC1 by Clinical Group")
ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, scatterpoints=1)
ax = axes[1]
pts = ax.scatter(scores["PC1"], scores["PC2"], c=metadata["ga_delivery_weeks"], cmap="Blues", s=70,
                 edgecolor=MUTED, linewidth=0.6, vmin=24, vmax=40)
fig.colorbar(pts, ax=ax, label="Weeks of pregnancy at delivery")
ax.set_xlabel(f"PC1 ({var[0]:.0%} of variation)"); ax.set_ylabel(f"PC2 ({var[1]:.0%} of variation)")
ax.set_title("Same Placentas Shaded by Delivery Week")
fig.tight_layout(); save(fig, "fig2_pca.png")

# %%
fig, ax = plt.subplots(figsize=(7.5, 4.4))
mlog = -np.log10(assoc.astype(float))
im = ax.imshow(mlog.values, cmap="Blues", vmin=0, vmax=5, aspect="auto")
ax.set_xticks(range(5)); ax.set_xticklabels([f"PC{i+1}\n({v:.0%})" for i, v in enumerate(var)])
ax.set_yticks(range(len(mlog))); ax.set_yticklabels(mlog.index)
for i in range(mlog.shape[0]):
    for j in range(mlog.shape[1]):
        p = assoc.iloc[i, j]
        ax.text(j, i, f"{p:.0e}" if p < 0.01 else f"{p:.2f}", ha="center", va="center", fontsize=8,
                color="white" if mlog.iloc[i, j] > 2.8 else INK)
fig.colorbar(im, ax=ax, label="Strength of link (-log10 p)")
ax.set_title("Which Clinical Factors Line Up With Each Main Pattern?")
ax.set_xlabel("Principal component (share of variation)")
fig.tight_layout(); save(fig, "fig3_pc_audit.png")

# %% [markdown]
# **Reading.** PC1, the largest pattern (22%), lines up with PE and not with delivery week; term and preterm controls sit
# together on the left. Delivery week lines up with PC2, together with growth restriction and sequencing depth. Fetal sex
# sits on PC5. First answer to Q6: the dominant PE pattern is not a timing pattern. The rest of Part B puts a number on it.

# %% [markdown]
# ## 4. Q6, Model Comparison on All 30 Placentas
#
# Four DESeq2 models, same thresholds, all contrasting PE vs no PE:
#
# | Model | Design | Role |
# |---|---|---|
# | Naive | `~ pe_status` | What an unadjusted comparison would report |
# | No timing | `~ fetal_sex + fgr_status + pe_status` | Reference for the timing question |
# | + weeks | `~ fetal_sex + fgr_status + gestational_age + pe_status` | Timing as a straight line in weeks |
# | + term indicator | `~ fetal_sex + fgr_status + term_birth + pe_status` | Timing as "term vs preterm", any shape |

# %%
def run_deseq(samples, design, contrast, meta=None, cnt=None, min_count=10, min_samples=5):
    meta = (metadata if meta is None else meta).loc[samples].copy()
    cnt = (counts if cnt is None else cnt)[samples]
    cnt = cnt[(cnt >= min_count).sum(axis=1) >= min(min_samples, len(samples))]
    if "ga_delivery_weeks" in meta:
        meta["ga_c"] = meta["ga_delivery_weeks"] - meta["ga_delivery_weeks"].mean()   # centred weeks
    for col in meta.columns:
        if meta[col].dtype == object:
            meta[col] = meta[col].astype(str)
    dds = DeseqDataSet(counts=cnt.T.astype(int), metadata=meta, design=design, quiet=True, inference=INF)
    dds.deseq2()
    st = DeseqStats(dds, contrast=contrast, quiet=True, inference=INF)
    st.summary()
    return st.results_df.sort_values("padj")

PADJ_CUTOFF, LFC_CUTOFF = 0.05, 1.0
def hits_of(res, padj=PADJ_CUTOFF, lfc=LFC_CUTOFF):
    return res[(res.padj < padj) & (res.log2FoldChange.abs() >= lfc)]

metadata["term_birth"] = np.where(metadata.group == "TB", "term", "preterm")
ALL = metadata.index
PE_C = ["pe_status", "PE", "no_PE"]
m_naive = run_deseq(ALL, "~pe_status", PE_C)
m_notime = run_deseq(ALL, "~fetal_sex + fgr_status + pe_status", PE_C)
m_weeks = run_deseq(ALL, "~fetal_sex + fgr_status + ga_c + pe_status", PE_C)
m_term = run_deseq(ALL, "~fetal_sex + fgr_status + term_birth + pe_status", PE_C)
for name, r in [("naive", m_naive), ("no_timing", m_notime), ("plus_weeks", m_weeks), ("plus_term_indicator", m_term)]:
    r.to_csv(RESULTS_DIR / f"de_all30_{name}.csv")
model_counts = pd.DataFrame({"model": ["Naive", "No timing (sex + FGR)", "+ weeks at delivery", "+ term indicator"],
                             "passing_genes": [len(hits_of(r)) for r in (m_naive, m_notime, m_weeks, m_term)]})
model_counts.to_csv(RESULTS_DIR / "q6_model_counts.csv", index=False); display(model_counts)

# %%
# Share of each gene's PE difference explained by delivery timing, for genes passing the no-timing model
ref = hits_of(m_notime).index
share = pd.DataFrame({
    "lfc_no_timing": m_notime.loc[ref, "log2FoldChange"],
    "lfc_plus_weeks": m_weeks.reindex(ref).log2FoldChange,
    "lfc_plus_term": m_term.reindex(ref).log2FoldChange})
share["share_weeks"] = 1 - share.lfc_plus_weeks / share.lfc_no_timing
share["share_term"] = 1 - share.lfc_plus_term / share.lfc_no_timing
share["still_pass_term"] = share.index.isin(hits_of(m_term).index)
share.to_csv(RESULTS_DIR / "q6_share_explained_by_timing.csv")
q6 = {}
for col, lab in [("share_weeks", "weeks"), ("share_term", "term indicator")]:
    s = share[col]
    q6[lab] = dict(median=s.median(), q25=s.quantile(.25), q75=s.quantile(.75), over_half=int((s > 0.5).sum()),
                   negative=int((s < 0).sum()), same_sign=(np.sign(share.lfc_no_timing) == np.sign(1 - s)).mean())
    print(f"Timing as {lab}: median share explained {s.median():.1%} (IQR {s.quantile(.25):.1%} to {s.quantile(.75):.1%}); "
          f"{(s > 0.5).sum()} of {len(s)} genes lose more than half their difference; {(s < 0).sum()} get larger.")
print(f"{share.still_pass_term.sum()} of {len(share)} genes still pass with the term indicator. "
      f"Direction kept for {(np.sign(share.lfc_plus_term) == np.sign(share.lfc_no_timing)).mean():.0%}.")
print("Unadjusted naive model: fetal sex + FGR adjustment alone changes the passing count from "
      f"{len(hits_of(m_naive))} to {len(hits_of(m_notime))}.")

# %%
KEY = ["FLT1", "ENG", "LEP", "HTRA4", "PAPPA2", "FSTL3", "INHBA", "SERPINA3"]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axes[0]
lim = [min(share.lfc_no_timing.min(), share.lfc_plus_term.min()) - 0.3, max(share.lfc_no_timing.max(), share.lfc_plus_term.max()) + 0.3]
ax.plot(lim, lim, color=MUTED, lw=1, ls="--", label="No change after adjustment")
ax.plot(lim, [0.5 * v for v in lim], color=GREY, lw=1, ls=":", label="Half explained by timing")
ax.scatter(share.lfc_no_timing, share.lfc_plus_term, s=18, color=C_PE, edgecolor="white", linewidth=0.4, zorder=3)
for g in ["FLT1", "ENG", "LEP", "HTRA4"]:
    if g in share.index:
        ax.annotate(g, (share.loc[g, "lfc_no_timing"], share.loc[g, "lfc_plus_term"]), xytext=(5, -10),
                    textcoords="offset points", fontsize=9)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("Preeclampsia fold change before accounting for timing (log2)")
ax.set_ylabel("Fold change after accounting for term vs preterm (log2)")
ax.set_title("PE Effects Before and After Timing Adjustment")
ax.legend(fontsize=8.5, loc="upper left")
ax = axes[1]
bins = np.linspace(-0.5, 0.8, 27)
ax.hist(share.share_term.clip(-0.5, 0.8) * 100, bins=bins * 100, color=C_PRE, alpha=0.85, label="Term vs preterm indicator")
ax.hist(share.share_weeks.clip(-0.5, 0.8) * 100, bins=bins * 100, histtype="step", lw=2, color=INK, label="Weeks at delivery")
ax.axvline(share.share_term.median() * 100, color=C_PRE, lw=1.5, ls="--")
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_xlabel("Change in estimated PE effect after timing adjustment (%)")
ax.set_ylabel("Number of genes")
ax.set_title(f"Median Share: {share.share_term.median():.0%} (Term) or {share.share_weeks.median():.0%} (Weeks)")
ax.legend(fontsize=8.5, loc="upper right")
fig.tight_layout(); save(fig, "fig4_timing_share.png")

# %% [markdown]
# ## 5. Q6, Variance Partition: How Much of Each Gene's Variation Tracks Timing vs Disease?
#
# For every expressed gene we fit a linear model of log2 CPM on PE, growth restriction, fetal sex and weeks at delivery
# (all 30 placentas) and record the extra variance each factor explains beyond the others (partial R squared).

# %%
X_df = pd.DataFrame({"pe": (metadata.pe_status == "PE").astype(float), "fgr": (metadata.fgr_status == "FGR").astype(float),
                     "male": (metadata.fetal_sex == "male").astype(float),
                     "ga": metadata.ga_delivery_weeks - metadata.ga_delivery_weeks.mean()}, index=metadata.index)
Y = log_cpm.loc[:, metadata.index].T.values
def rss(cols):
    X = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in cols])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return ((Y - X @ beta) ** 2).sum(axis=0)
tss = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
full = rss(["pe", "fgr", "male", "ga"])
vp = pd.DataFrame({k: (rss([c for c in ["pe", "fgr", "male", "ga"] if c != k]) - full) / tss for k in ["pe", "fgr", "male", "ga"]},
                  index=log_cpm.index)
vp.columns = ["Preeclampsia", "Growth restriction", "Fetal sex", "Weeks at delivery"]
vp.to_csv(RESULTS_DIR / "q6_variance_partition.csv")
vp_sets = {"All expressed genes": vp, "Genes passing (no-timing model)": vp.loc[vp.index.intersection(ref)]}
vp_summary = pd.DataFrame({k: v.median() for k, v in vp_sets.items()}).T
display((vp_summary * 100).round(1))

# %%
fig, ax = plt.subplots(figsize=(9, 4.6))
d = vp_sets["Genes passing (no-timing model)"] * 100
parts = ax.boxplot([d[c] for c in d.columns], vert=False, widths=0.55, patch_artist=True, showfliers=False,
                   medianprops=dict(color=INK, lw=1.5))
for p, col in zip(parts["boxes"], [C_PE, C_PRE, GREY, C_TERM]):
    p.set_facecolor(col); p.set_alpha(0.85); p.set_edgecolor("white")
ax.set_yticks(range(1, 5)); ax.set_yticklabels(d.columns)
ax.set_xlabel("Extra share of each gene's variation explained, beyond the other factors (%)")
ax.set_title(f"Partial R-Squared Across {len(d)} Selected Genes")
fig.tight_layout(); save(fig, "fig5_variance_partition.png")

# %% [markdown]
# ## 6. Q6, Preterm-Only Adjusted Comparison
#
# This secondary model compares 15 PE placentas with 10 non-PE preterm placentas (5 spontaneous preterm and 5 FGR-only),
# adjusting for weeks, FGR and fetal sex. The groups overlap in delivery weeks but are not individually matched.
# It supplies the gene list for Parts C to E. The 5 term vs 5 spontaneous-preterm control comparison is a timing-related
# contrast, not a pure timing effect: delivery circumstances and other differences may also contribute.

# %%
preterm = metadata.index[metadata.group != "TB"]
print(metadata.loc[preterm].groupby(["pe_status", "fgr_status"]).size().to_string())
primary = run_deseq(preterm, "~fetal_sex + ga_c + fgr_status + pe_status", PE_C)
primary.to_csv(RESULTS_DIR / "de_matched_PE_vs_noPE_preterm.csv")
hits = hits_of(primary); hits.to_csv(RESULTS_DIR / "passing_genes_matched.csv")
up, down = hits.index[hits.log2FoldChange > 0], hits.index[hits.log2FoldChange < 0]
print(f"\nTiming-matched model: tested {len(primary):,} genes; {len(hits)} pass ({len(up)} higher in PE, {len(down)} lower).")

controls = metadata.index[metadata.group.isin(["TB", "iPTB"])]
ga_sig = run_deseq(controls, "~fetal_sex + term_birth", ["term_birth", "term", "preterm"])
ga_sig.to_csv(RESULTS_DIR / "de_timing_signature_term_vs_preterm_controls.csv")
h_ga = hits_of(ga_sig)
common = primary.index.intersection(m_notime.index)
r_match = stats.pearsonr(primary.loc[common, "log2FoldChange"], m_notime.loc[common, "log2FoldChange"])[0]
print(f"Timing signature (term vs preterm controls): {len(h_ga)} genes.")
print(f"Overlap of timing signature with timing-matched PE genes: {len(set(h_ga.index) & set(hits.index))} genes.")
print(f"Fold changes, timing-matched model vs all-30 no-timing model, all tested genes: r = {r_match:.2f}")
in_both = hits.index.intersection(hits_of(m_notime).index)
print(f"{len(in_both)} of {len(hits)} timing-matched genes also pass in the all-30 no-timing model.")

# %% [markdown]
# ### 6b. Is the small timing signature just a small-sample effect? And what would a term-only study conclude?
#
# Two checks contributed by a teammate's parallel analysis:
#
# * **Noise and power baseline.** The timing signature compares 5 vs 5 placentas. We repeat 5-vs-5 tests 8 times each:
#   (i) the 10 control placentas split at random (no real difference, so this measures noise), and (ii) 5 randomly chosen PE
#   placentas vs the 5 spontaneous preterm births (a real disease difference at the same sample size). Same design and
#   thresholds as the timing signature.
# * **Worst case.** Many studies only have term controls. One model with three groups (PE, spontaneous preterm, term) splits
#   each gene's PE-vs-term difference exactly into a disease part (PE vs preterm) and a timing part (preterm vs term). For
#   genes passing PE vs term, we report how much of the difference is the timing part.

# %%
def five_v_five(a_ids, b_ids):
    m = metadata.loc[list(a_ids) + list(b_ids), ["fetal_sex"]].copy()
    m["grp"] = ["A"] * len(a_ids) + ["B"] * len(b_ids)
    return len(hits_of(run_deseq(m.index, "~fetal_sex + grp", ["grp", "A", "B"], meta=m)))
iptb_ids = metadata.index[metadata.group == "iPTB"].tolist(); tb_ids = metadata.index[metadata.group == "TB"].tolist()
pe_ids = metadata.index[metadata.pe_status == "PE"].tolist()
noise_draws, disease_draws = [], []
for _ in range(8):
    split = list(rng.permutation(iptb_ids + tb_ids))
    while set(split[:5]) in (set(iptb_ids), set(tb_ids)):
        split = list(rng.permutation(iptb_ids + tb_ids))
    noise_draws.append(five_v_five(split[:5], split[5:]))
    disease_draws.append(five_v_five(rng.choice(pe_ids, 5, replace=False), iptb_ids))
power = pd.DataFrame({"random_split_of_controls": noise_draws, "five_PE_vs_five_preterm": disease_draws})
power.to_csv(RESULTS_DIR / "q6_noise_power_check.csv", index=False)
print("Random splits of controls (noise):  ", sorted(noise_draws))
print("5 PE vs 5 spontaneous preterm:      ", sorted(disease_draws))
print("Term vs preterm controls (timing):  ", len(h_ga))

# %%
three = metadata.index[metadata.group.isin(["EOPE", "EOPE+FGR", "iPTB", "TB"])]
m3 = metadata.loc[three, ["fetal_sex", "fgr_status"]].copy()
m3["grp3"] = np.select([metadata.loc[three, "pe_status"] == "PE", metadata.loc[three, "group"] == "iPTB"], ["PE", "preterm"], "term")
dds3 = DeseqDataSet(counts=counts.loc[expressed, three].T.astype(int), metadata=m3, design="~fetal_sex + fgr_status + grp3",
                    quiet=True, inference=INF)
dds3.deseq2()
def contrast3(a, b):
    st = DeseqStats(dds3, contrast=["grp3", a, b], quiet=True, inference=INF); st.summary(); return st.results_df
pe_term, pe_pre, pre_term = contrast3("PE", "term"), contrast3("PE", "preterm"), contrast3("preterm", "term")
wc_hits = hits_of(pe_term).index
timing_part = (pre_term.loc[wc_hits, "log2FoldChange"] / pe_term.loc[wc_hits, "log2FoldChange"])
holds = (pe_pre.loc[wc_hits, "padj"] < 0.05) & (np.sign(pe_pre.loc[wc_hits, "log2FoldChange"]) == np.sign(pe_term.loc[wc_hits, "log2FoldChange"]))
worst = pd.DataFrame({"pe_vs_term": pe_term.loc[wc_hits, "log2FoldChange"], "pe_vs_preterm": pe_pre.loc[wc_hits, "log2FoldChange"],
                      "preterm_vs_term": pre_term.loc[wc_hits, "log2FoldChange"], "timing_part": timing_part, "holds_vs_preterm": holds})
worst.to_csv(RESULTS_DIR / "q6_worst_case_term_only.csv")
print(f"Term-only design: {len(wc_hits)} genes pass PE vs term. Median timing part = {timing_part.median():.0%} of the difference. "
      f"{holds.mean():.0%} still pass against spontaneous preterm controls; "
      f"{((~holds) & (timing_part >= 0.5)).mean():.0%} look mostly timing-driven.")
worst_case_median = timing_part.median()

# %%
grid = pd.DataFrame({f"|LFC|>={l}": [len(hits_of(primary, p, l)) for p in (0.01, 0.05, 0.1)] for l in (0, 0.58, 1, 2)},
                    index=["padj<0.01", "padj<0.05", "padj<0.10"])
grid.to_csv(RESULTS_DIR / "threshold_sensitivity.csv"); display(grid)

# %%
LABEL = ["FLT1", "ENG", "LEP", "HTRA4", "PAPPA2", "FSTL3", "SERPINA3", "INHBA", "NDRG1", "HSD17B1"]
fig, ax = plt.subplots(figsize=(8, 5.6))
y = -np.log10(primary.padj.clip(lower=1e-300))
ax.scatter(primary.log2FoldChange, y, s=5, color=GREY, rasterized=True, label="Not passing")
ax.scatter(primary.loc[up, "log2FoldChange"], y[up], s=12, color=C_PE, label=f"Higher in preeclampsia ({len(up)})")
ax.scatter(primary.loc[down, "log2FoldChange"], y[down], s=12, color=C_PRE, label=f"Lower in preeclampsia ({len(down)})")
ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", lw=0.8, color=MUTED)
for x0 in (-LFC_CUTOFF, LFC_CUTOFF):
    ax.axvline(x0, ls="--", lw=0.8, color=MUTED)
for g in LABEL:
    if g in hits.index:
        ax.annotate(g, (primary.loc[g, "log2FoldChange"], y[g]), xytext=(4, 2), textcoords="offset points", fontsize=8.5)
ax.set_xlabel("Fold change, preeclampsia vs. no preeclampsia (log2 scale; +1 = twice as high)")
ax.set_ylabel("Strength of evidence (-log10 adjusted p)")
ax.set_title(f"Adjusted Preterm Comparison: {len(hits)} Genes Pass")
ax.legend(fontsize=9, loc="upper left")
fig.tight_layout(); save(fig, "fig6_volcano.png")

# %%
# Key genes: fold change with 95% Wald confidence interval
key = primary.reindex(KEY + ["HSD17B1"])[["baseMean", "log2FoldChange", "lfcSE", "padj"]].copy()
key["fold"] = 2 ** key.log2FoldChange
key["fold_low"] = 2 ** (key.log2FoldChange - 1.96 * key.lfcSE); key["fold_high"] = 2 ** (key.log2FoldChange + 1.96 * key.lfcSE)
key["share_explained_by_timing_term"] = share.reindex(key.index).share_term
key.to_csv(RESULTS_DIR / "key_genes.csv"); display(key.round(3))

# %% [markdown]
# ### Pathway context (Hallmark gene sets)
#
# Fisher's exact test of passing genes against all tested genes for the 50 MSigDB Hallmark sets, BH-adjusted. MSigDB's
# own download now requires sign-in, so the notebook fetches the Hallmark v7.0 file from the GSEApy test-data mirror and
# verifies its checksum. A local `data/track01_gene_sets.gmt` is accepted only if it matches the same v7.0 checksum;
# the official `get_gene_sets.py` downloads a different release and cannot be substituted without changing the analysis.

# %%
GMT_URL = "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt"
GMT_SHA = "4803c09aae50bbd96d53fc681bd5545a4c89278993c15145dd7b7c994f893516"
gmt_path = DATA_DIR / "track01_gene_sets.gmt"
if not gmt_path.exists():
    gmt_path = RESULTS_DIR / "h.all.v7.0.symbols.gmt"
    if not gmt_path.exists():
        from urllib.request import urlopen
        gmt_path.write_bytes(urlopen(GMT_URL).read())
if hashlib.sha256(gmt_path.read_bytes()).hexdigest() != GMT_SHA:
    raise ValueError(f"{gmt_path} is not the expected Hallmark v7.0 file (SHA-256 {GMT_SHA}). "
                     "Use the matching v7.0 file; the 2024.1 release from get_gene_sets.py is different.")
hallmark = {ln.split("\t")[0].replace("HALLMARK_", ""): set(ln.strip().split("\t")[2:])
            for ln in gmt_path.read_text().splitlines() if ln.strip()}

def ora(genes, universe):
    genes, universe = set(genes), set(universe); rows = []
    for name, gs in hallmark.items():
        gs = gs & universe
        a = len(genes & gs); b = len(genes) - a; c_ = len(gs) - a; d = len(universe) - a - b - c_
        orr, p = stats.fisher_exact([[a, b], [c_, d]], alternative="greater")
        rows.append((name, a, len(gs), orr, p, ", ".join(sorted(genes & gs))))
    out = pd.DataFrame(rows, columns=["set", "overlap", "set_size", "odds_ratio", "p", "genes"])
    out["padj"] = multipletests(out.p, method="fdr_bh")[1]
    return out.sort_values("p")
ora_up = ora(up, primary.index); ora_down = ora(down, primary.index); ora_timing = ora(h_ga.index, ga_sig.index)
ora_up.to_csv(RESULTS_DIR / "hallmark_up.csv", index=False); ora_down.to_csv(RESULTS_DIR / "hallmark_down.csv", index=False)
ora_timing.to_csv(RESULTS_DIR / "hallmark_timing_signature.csv", index=False)
display(ora_up.head(6)); display(ora_down.head(3)); display(ora_timing.head(3))

# %%
top = ora_up.head(6).iloc[::-1]
fig, ax = plt.subplots(figsize=(8, 4))
ax.barh(top.set.str.replace("_", " ").str.title().str.replace("Mtorc1", "mTORC1"), -np.log10(top.padj), color=C_PE, height=0.6)
for i, (o, s) in enumerate(zip(top.overlap, top.set_size)):
    ax.text(-np.log10(top.padj.iloc[i]) + 0.05, i, f"{o} of {s} genes", va="center", fontsize=9, color=MUTED)
ax.axvline(-np.log10(0.05), ls="--", lw=0.8, color=MUTED)
ax.set_xlabel("Enrichment among genes higher in preeclampsia (-log10 adjusted p)")
ax.set_title("Pathway Enrichment Among Genes Higher in PE")
fig.tight_layout(); save(fig, "fig7_pathways.png")

# %% [markdown]
# ## 7. Robustness of the Timing-Independent Signal
#
# * **(a)** drop the samples with blood-pressure audit flags or low read assignment;
# * **(b)** controls = spontaneous preterm births only;
# * **(c)** add maternal age and sequencing depth;
# * **(d)** controls = growth restriction only. Spontaneous preterm births are usually laboured vaginal deliveries, while PE
#   and FGR deliveries are mostly pre-labour caesareans. Labour is not recorded in the bulk data, so comparing PE with FGR
#   placentas is our best available proxy for matching on labour and delivery mode.
# * Leave-one-out over the 10 preterm non-PE controls, and 10 label permutations within FGR strata.

# %%
DESIGN = "~fetal_sex + ga_c + fgr_status + pe_status"
drop = sorted(set(flagged_bp) | set(low_qc))
print("Sensitivity (a) drops:", drop)
sens = {
    "(a) drop flagged / low-QC": run_deseq(preterm.difference(drop), DESIGN, PE_C),
    "(b) controls = spontaneous preterm only": run_deseq(metadata.index[metadata.group.isin(["EOPE", "EOPE+FGR", "iPTB"])], DESIGN, PE_C),
}
meta_c = metadata.copy(); meta_c["age_c"] = meta_c.maternal_age_years - meta_c.maternal_age_years.mean()
meta_c["depth_c"] = np.log2(meta_c.library_M) - np.log2(meta_c.library_M).mean()
sens["(c) + maternal age + depth"] = run_deseq(preterm, "~fetal_sex + ga_c + age_c + depth_c + fgr_status + pe_status", PE_C, meta=meta_c)
sens["(d) controls = growth restriction only"] = run_deseq(metadata.index[metadata.group.isin(["EOPE", "EOPE+FGR", "FGR"])], DESIGN, PE_C)
rows = []
for k, r in sens.items():
    r.to_csv(RESULTS_DIR / f"de_sensitivity_{k[1]}.csv")
    cm_ = hits.index.intersection(r.index)
    rows.append(dict(model=k, passing=len(hits_of(r)),
                     same_direction=(np.sign(r.loc[cm_, "log2FoldChange"]) == np.sign(hits.loc[cm_, "log2FoldChange"])).mean(),
                     still_padj_lt_0_05=(r.loc[cm_, "padj"] < 0.05).mean(),
                     lfc_correlation=stats.pearsonr(r.loc[cm_, "log2FoldChange"], hits.loc[cm_, "log2FoldChange"])[0],
                     median_effect_ratio=(r.loc[cm_, "log2FoldChange"] / hits.loc[cm_, "log2FoldChange"]).median()))
sens_tab = pd.DataFrame(rows); sens_tab.to_csv(RESULTS_DIR / "sensitivity_models.csv", index=False)
display(sens_tab.round(2))

# %%
loo_hits, perm_counts = {}, []
if RUN_SLOW:
    pre_controls = metadata.index[(metadata.group != "TB") & (metadata.pe_status == "no_PE")]
    for s in pre_controls:
        loo_hits[s] = set(hits_of(run_deseq(preterm.drop(s), DESIGN, PE_C)).index)
        print(f"  without {s}: {len(loo_hits[s])} pass")
    for i in range(10):
        m_perm = metadata.loc[preterm].copy()
        for fgr in ("FGR", "no_FGR"):
            idx = m_perm.index[m_perm.fgr_status == fgr]
            m_perm.loc[idx, "pe_status"] = rng.permutation(m_perm.loc[idx, "pe_status"].values)
        perm_counts.append(len(hits_of(run_deseq(preterm, DESIGN, PE_C, meta=m_perm))))
    print("Permutation pass counts:", perm_counts)
    retention = pd.Series({g: np.mean([g in v for v in loo_hits.values()]) for g in hits.index}, name="loo_retention")
    retention.to_csv(RESULTS_DIR / "loo_retention.csv")
    print(f"{(retention == 1).sum()} of {len(hits)} genes pass in all 10 leave-one-out refits; {(retention >= 0.8).sum()} in at least 8.")
else:
    retention = pd.Series(np.nan, index=hits.index, name="loo_retention")

# %%
if RUN_SLOW:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    ax.bar(range(len(loo_hits)), [len(v) for v in loo_hits.values()], color=C_PRE, width=0.6)
    ax.axhline(len(hits), color=C_PE, lw=2, label=f"All samples ({len(hits)})")
    ax.set_xticks(range(len(loo_hits))); ax.set_xticklabels([s.replace("_VT", "") for s in loo_hits], rotation=45)
    ax.set_xlabel("Control pregnancy left out"); ax.set_ylabel("Genes passing")
    ax.set_title("Leave-One-Control-Out Sensitivity"); ax.legend(fontsize=9, loc="lower right")
    ax = axes[1]
    ax.bar(["Real labels"] + [f"Shuffle {i+1}" for i in range(len(perm_counts))], [len(hits)] + perm_counts,
           color=[C_PE] + [GREY] * len(perm_counts), width=0.6)
    ax.set_ylabel("Genes passing"); ax.tick_params(axis="x", rotation=60)
    ax.set_title("Shuffled Labels Give Almost No Genes")
    fig.tight_layout(); save(fig, "fig8_robustness.png")

# %% [markdown]
# ## 8. Q9, Cell Mix in the Bulk Tissue (Marker Scores)
#
# For each broad class in the marker table: mean z-scored log CPM of its markers per placenta, then PE vs no PE with the
# timing-matched covariates (linear regression, preterm placentas, 95% CI). If a class's markers rise together, a change
# in that cell type's share of the tissue would be the simplest explanation for its genes rising.

# %%
z = log_cpm_all.sub(log_cpm_all.mean(axis=1), axis=0).div(log_cpm_all.std(axis=1) + 1e-9, axis=0)
present = markers[markers.gene_symbol.isin(counts.index)]
class_scores = pd.DataFrame({cls: z.loc[g.gene_symbol].mean() for cls, g in present.groupby("cell_class")})
dfm = metadata.loc[preterm, ["pe_status", "fgr_status", "fetal_sex", "ga_delivery_weeks"]].join(class_scores.loc[preterm])
mrows = []
for cls in class_scores.columns:
    fit = smf.ols(f"{cls} ~ C(pe_status, Treatment('no_PE')) + C(fgr_status) + C(fetal_sex) + ga_delivery_weeks", dfm).fit()
    t = [x for x in fit.params.index if x.startswith("C(pe_status")][0]; ci = fit.conf_int().loc[t]
    mrows.append(dict(cell_class=cls, n_markers=int((present.cell_class == cls).sum()), pe_effect_sd=fit.params[t],
                      ci_low=ci[0], ci_high=ci[1], p=fit.pvalues[t]))
mix = pd.DataFrame(mrows); mix["padj"] = multipletests(mix.p, method="fdr_bh")[1]
mix = mix.sort_values("pe_effect_sd"); mix.to_csv(RESULTS_DIR / "marker_class_scores.csv", index=False)
display(mix.round(3))

# %% [markdown]
# ## 9. Q9, Are the Genes Normally Broad or Concentrated? (Healthy Single-Cell Reference)
#
# Admati/Skarbianskis pseudobulk, **control donors only** (healthy reference, 9 donors), donor x class columns with at
# least 10 cells. For each gene: mean CPM per class across donors, then the tau specificity index on log2(CPM + 1).
# Genes expressed (> 5 CPM) in at least one class are kept. The null compares each gene set with 1,000 random sets drawn
# to match its bulk expression level (deciles of baseMean).

# %%
PB = DATA_DIR / "extension" / "pseudobulk"
pb_counts_raw = pd.read_csv(PB / "track01_pseudobulk_counts.csv", index_col=0)
pb_meta_all = pd.read_csv(PB / "track01_pseudobulk_metadata.csv", index_col=0)
donors = pd.read_csv(PB / "track01_pseudobulk_donors.csv", index_col=0)
pb_meta = pb_meta_all[pb_meta_all.n_cells >= 10]
pb_counts = pb_counts_raw[pb_meta.index]
pb_cpm = pb_counts.div(pb_counts.sum(axis=0), axis=1) * 1e6
pb_log = np.log2(pb_cpm + 1)
print(f"Pseudobulk: {pb_counts.shape[0]:,} genes x {pb_counts.shape[1]} donor-by-class columns from {pb_meta.donor_id.nunique()} donors.")
ctrl_cols = pb_meta.index[pb_meta.pe_status == "control"]
print("Healthy reference: %d control donors; columns per class:" % pb_meta.loc[ctrl_cols, "donor_id"].nunique())
print(pb_meta.loc[ctrl_cols].cell_class.value_counts().to_string())
healthy = np.log2(pb_cpm[ctrl_cols].T.groupby(pb_meta.loc[ctrl_cols, "cell_class"]).mean().T + 1)
healthy = healthy[healthy.max(axis=1) > np.log2(6)]
tau = (1 - healthy.div(healthy.max(axis=1), axis=0)).sum(axis=1) / (healthy.shape[1] - 1)
top_class = healthy.idxmax(axis=1)

bg = primary.loc[primary.index.intersection(tau.index), ["baseMean"]].copy()
bg["bin"] = pd.qcut(np.log10(bg.baseMean + 1), 10, labels=False)
def matched_null(genes, n=1000):
    genes = [g for g in genes if g in bg.index]; need = bg.loc[genes, "bin"].value_counts()
    return np.array([np.median(tau[np.concatenate([rng.choice(bg.index[bg.bin == b], k, replace=False) for b, k in need.items()])])
                     for _ in range(n)])
timing_genes = [g for g in h_ga.index if g in tau.index]
spec_rows = []
for name, gs in [("Higher in PE", up), ("Lower in PE", down), ("Delivery-timing signature", timing_genes)]:
    gs = [g for g in gs if g in tau.index]
    nul = matched_null(gs)
    spec_rows.append(dict(gene_set=name, n=len(gs), median_tau=tau[gs].median(), null_median_tau=np.median(nul),
                          empirical_p=(np.sum(nul >= tau[gs].median()) + 1) / (len(nul) + 1),
                          share_highly_specific=(tau[gs] > 0.8).mean()))
spec_rows.append(dict(gene_set="All expressed genes", n=len(bg), median_tau=tau[bg.index].median(), null_median_tau=np.nan,
                      empirical_p=np.nan, share_highly_specific=(tau[bg.index] > 0.8).mean()))
spec = pd.DataFrame(spec_rows); spec.to_csv(RESULTS_DIR / "q9_specificity.csv", index=False); display(spec.round(3))

loc_tab = pd.DataFrame({"higher_in_PE": top_class.reindex([g for g in up if g in tau.index]).value_counts(),
                        "lower_in_PE": top_class.reindex([g for g in down if g in tau.index]).value_counts(),
                        "timing_signature": top_class.reindex(timing_genes).value_counts(),
                        "all_expressed": top_class.reindex(bg.index).value_counts()}).fillna(0).astype(int)
loc_tab.to_csv(RESULTS_DIR / "q9_top_cell_class.csv"); display(loc_tab.sort_values("higher_in_PE", ascending=False))
troph = ["syncytiotrophoblast", "extravillous_trophoblast", "villous_cytotrophoblast"]
n_up_t = int(loc_tab.loc[loc_tab.index.isin(troph), "higher_in_PE"].sum()); n_up = int(loc_tab.higher_in_PE.sum())
bg_t = loc_tab.loc[loc_tab.index.isin(troph), "all_expressed"].sum() / loc_tab.all_expressed.sum()
p_troph = stats.binomtest(n_up_t, n_up, bg_t, alternative="greater").pvalue
print(f"{n_up_t} of {n_up} genes higher in PE peak in a trophoblast class in healthy placenta, vs {bg_t:.0%} of all "
      f"expressed genes (binomial p = {p_troph:.1e}).")
timing_t = int(loc_tab.loc[loc_tab.index.isin(troph), "timing_signature"].sum()); timing_n = int(loc_tab.timing_signature.sum())
print(f"Delivery-timing genes: {timing_t} of {timing_n} peak in trophoblast (binomial p vs background = "
      f"{stats.binomtest(timing_t, timing_n, bg_t, alternative='greater').pvalue:.2f}).")
n_ctrl_evt = pb_meta.loc[ctrl_cols].query("cell_class=='extravillous_trophoblast'").shape[0]
print(f"Note: the healthy extravillous trophoblast profile rests on {n_ctrl_evt} control donors.")

# %%
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axes[0]
sets = [("All expressed genes", list(bg.index), GREY), ("Delivery-timing\nsignature", timing_genes, C_TERM),
        ("Lower in\npreeclampsia", [g for g in down if g in tau.index], C_PRE), ("Higher in\npreeclampsia", [g for g in up if g in tau.index], C_PE)]
bp = ax.boxplot([tau[g] for _, g, _ in sets], vert=False, widths=0.55, patch_artist=True, showfliers=False,
                medianprops=dict(color=INK, lw=1.5))
for b, (_, _, col) in zip(bp["boxes"], sets):
    b.set_facecolor(col); b.set_edgecolor("white")
for i, (_, g, col) in enumerate(sets[1:], start=2):
    ax.scatter(tau[g], i + rng.uniform(-0.18, 0.18, len(g)), s=8, color=INK, alpha=0.35, zorder=3)
ax.set_yticks(range(1, 5)); ax.set_yticklabels([s[0] + f"\n(n={len(s[1])})" for s in sets], fontsize=9)
ax.set_xlabel("Cell-type specificity in healthy placenta\n(0 = made equally by all cell types, 1 = only one)")
ax.set_title("Preeclampsia Genes Are Concentrated, Not Broad")
ax = axes[1]
show = loc_tab.sort_values("higher_in_PE", ascending=True).tail(7)
yy = np.arange(len(show))
ax.barh(yy + 0.2, show.higher_in_PE / loc_tab.higher_in_PE.sum() * 100,
        height=0.38, color=C_PE, label="Genes higher in preeclampsia")
ax.barh(yy - 0.2, show.all_expressed / loc_tab.all_expressed.sum() * 100, height=0.38, color=GREY, label="All expressed genes")
ax.set_yticks(yy); ax.set_yticklabels([nice(c) for c in show.index])
ax.set_xlabel("Share of genes whose highest cell type is this one (%)")
ax.set_title("Where They Peak: Trophoblast")
ax.legend(fontsize=9, loc="lower right")
fig.tight_layout(); save(fig, "fig9_specificity.png")

# %% [markdown]
# ## 10. Q9, In Which Cell Type Do the Genes Change With Preeclampsia?
#
# Within each class we compare PE and control donors per gene (log2 CPM, donor = unit) over the bulk genes higher in PE
# that the class expresses (> 5 CPM). Each class is compared with 500 random gene sets of the same size drawn from genes
# that class expresses. **Late onset** (7 PE vs 6 controls, all delivered 37 to 40 weeks) is matched on delivery timing by
# design; **early onset** (10 PE vs 3 controls) is exploratory as instructed.
#
# The single-cell donor table records delivery mode, labour induction and magnesium sulfate, which the bulk data lack.
# In the late arm we therefore repeat the syncytiotrophoblast comparison restricted to caesarean deliveries and to donors
# who did not receive magnesium sulfate.

# %%
print(pd.crosstab([donors.onset_group, donors.pe_status], donors.delivery_mode).to_string())
print(pd.crosstab([donors.onset_group, donors.pe_status], donors.magnesium_sulfate).to_string())

def within_class(cls, onset, donor_filter=None, n_draw=500):
    cols = pb_meta.index[(pb_meta.cell_class == cls) & (pb_meta.onset_group == onset)]
    if donor_filter is not None:
        cols = cols[pb_meta.loc[cols, "donor_id"].isin(donor_filter)]
    pe = cols[pb_meta.loc[cols, "pe_status"] == "PE"]; ct = cols[pb_meta.loc[cols, "pe_status"] == "control"]
    if len(pe) < 2 or len(ct) < 2:
        return None
    pool = pb_cpm.index[pb_cpm[cols].mean(axis=1) > 5]
    diff = pb_log.loc[pool, pe].mean(axis=1) - pb_log.loc[pool, ct].mean(axis=1)
    genes = [g for g in up if g in diff.index]
    obs = diff[genes].median()
    draws = np.array([np.median(rng.choice(diff.values, len(genes), replace=False)) for _ in range(n_draw)])
    return dict(cell_class=cls, onset=onset, n_pe=len(pe), n_ctrl=len(ct), genes_expressed=len(genes),
                frac_higher_in_pe=(diff[genes] > 0).mean(), median_log2_diff=obs, null_p=(np.sum(draws >= obs) + 1) / (n_draw + 1))
conc = pd.DataFrame([r for c in sorted(pb_meta.cell_class.unique()) for o in ("late", "early") if (r := within_class(c, o))])
conc.to_csv(RESULTS_DIR / "q9_within_class_concordance.csv", index=False)
display(conc.query("onset=='late'").sort_values("median_log2_diff", ascending=False).round(3))
display(conc.query("onset=='early'").sort_values("median_log2_diff", ascending=False).round(3))

late_csec = donors.index[(donors.onset_group == "late") & (donors.delivery_mode == "c_section")]
late_nomg = donors.index[(donors.onset_group == "late") & (donors.magnesium_sulfate == "no")]
stb_checks = pd.DataFrame([dict(subset="Late onset, all donors", **within_class("syncytiotrophoblast", "late")),
                           dict(subset="Late onset, caesarean only", **within_class("syncytiotrophoblast", "late", late_csec)),
                           dict(subset="Late onset, no magnesium sulfate", **within_class("syncytiotrophoblast", "late", late_nomg))])
stb_checks.to_csv(RESULTS_DIR / "q9_stb_delivery_checks.csv", index=False); display(stb_checks.round(3))

# %%
# Formal DESeq2 inside syncytiotrophoblast: late onset with delivery mode as a covariate, and early onset (exploratory)
pbm = pb_meta   # already carries delivery_mode for each donor
stb_late = pbm.index[(pbm.cell_class == "syncytiotrophoblast") & (pbm.onset_group == "late")]
stb_early = pbm.index[(pbm.cell_class == "syncytiotrophoblast") & (pbm.onset_group == "early")]
stb_de_late = run_deseq(stb_late, "~delivery_mode + pe_status", ["pe_status", "PE", "control"],
                        meta=pbm[["pe_status", "delivery_mode"]], cnt=pb_counts, min_count=5, min_samples=4)
stb_de_early = run_deseq(stb_early, "~pe_status", ["pe_status", "PE", "control"], meta=pbm[["pe_status"]],
                         cnt=pb_counts, min_count=5, min_samples=3)
stb_de_late.to_csv(RESULTS_DIR / "de_pseudobulk_STB_late_deliverymode.csv"); stb_de_early.to_csv(RESULTS_DIR / "de_pseudobulk_STB_early.csv")
for name, r in [("late (7 vs 6, + delivery mode)", stb_de_late), ("early (10 vs 3, exploratory)", stb_de_early)]:
    g = [x for x in up if x in r.index]
    print(f"STB {name}: {len(g)} bulk up-genes tested; {(r.loc[g,'log2FoldChange']>0).mean():.0%} higher in PE; "
          f"{(r.loc[g,'padj']<0.05).sum()} reach padj<0.05 on their own")

# %%
order = conc.query("onset=='late'").sort_values("median_log2_diff").cell_class.tolist()
order = [c for c in order if conc.query("onset=='late' and cell_class==@c").n_ctrl.iloc[0] >= 4]
fig, ax = plt.subplots(figsize=(8.5, 5.4))
for onset, col, mk, off in [("late", C_PRE, "s", 0.17), ("early", C_TERM, "^", -0.17)]:
    d = conc.set_index(["onset", "cell_class"]).loc[onset].reindex(order)
    ax.scatter(d.median_log2_diff, np.arange(len(order)) + off, color=col, marker=mk, s=60, zorder=3,
               label="Late onset: 7 PE vs 6 controls, same delivery weeks" if onset == "late" else "Early onset: 10 PE vs 3 controls (exploratory)")
ax.axvline(0, color=MUTED, lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels([nice(c) for c in order])
ax.set_xlabel("Median change in preeclampsia across our bulk 'higher' genes\n(log2, within each cell type, donor = unit)")
ax.set_title("The Genes Rise Most Inside Syncytiotrophoblast")
ax.legend(fontsize=8.5, loc="lower right", scatterpoints=1)
fig.tight_layout(); save(fig, "fig10_within_celltype.png")

# %%
show_classes = ["syncytiotrophoblast", "villous_cytotrophoblast", "extravillous_trophoblast", "stromal_cell",
                "endothelial_cell", "hofbauer_cell", "macrophage", "nk_cell", "t_cell"]
genes_plot = [g for g in KEY + ["HSD17B1"] if g in healthy.index]
heat = healthy.loc[genes_plot, show_classes]; heat = heat.div(heat.max(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8.5, 5))
im = ax.imshow(heat.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_yticks(range(len(genes_plot))); ax.set_yticklabels([f"{g} (tau {tau[g]:.2f})" for g in genes_plot])
ax.set_xticks(range(len(show_classes))); ax.set_xticklabels([nice(c) for c in show_classes], rotation=40, ha="right", fontsize=9)
fig.colorbar(im, ax=ax, label="Expression relative to the gene's\nhighest cell type (healthy donors)", shrink=0.8)
ax.set_title("Key Genes in Healthy Placenta: Mostly Trophoblast")
fig.tight_layout(); save(fig, "fig11_celltype_heatmap.png")

# %% [markdown]
# ## 10b. Gene Browser Lookups (Zeisel/Admati)
#
# Recorded on 14 November 2026 from the [Zeisel/Admati browser](https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html)
# for seven genes: five higher in PE, one lower (HSD17B1) and one timing gene (HLA-G). The browser marks PE vs control
# differences per cell, which treats cells as replicates, so we use it as qualitative context only. The same file is read
# by the workflow notebook.

# %%
browser_lookups = pd.read_csv(Path("inputs/track01_browser_lookups.csv"))
required = {"gene", "cell_types_expressing", "specificity", "differs_PE_vs_control", "notes"}
assert required <= set(browser_lookups.columns), f"missing columns: {required - set(browser_lookups.columns)}"
assert browser_lookups.gene.is_unique and len(browser_lookups) >= 3
browser_lookups.to_csv(RESULTS_DIR / "browser_lookups.csv", index=False)
display(browser_lookups)

# %% [markdown]
# ## 11. Headline: A Cross-Validated Low-Oxygen Score Against Delivery Week
#
# A per-placenta score for the hypoxia-response genes that rise in PE. To avoid circularity, the genes are **chosen
# without the placenta being scored**: for each of the 25 preterm placentas we refit a linear model on the other 24
# (log2 CPM ~ PE + FGR + sex + weeks), keep Hallmark hypoxia genes with PE effect >= 1 log2 and p < 0.001, and score the
# held-out placenta against the training non-PE mean. Term births were never in any training set and are scored with the
# genes chosen from all 25. If the PE gap were an artefact of gene selection, it would collapse here.

# %%
hyp_in = [g for g in hallmark["HYPOXIA"] if g in log_cpm.index]
Lh = log_cpm.loc[hyp_in]
def select_genes(train):
    Xt = np.column_stack([np.ones(len(train)), X_df.loc[train, ["pe", "fgr", "male", "ga"]].values])
    Yt = Lh[train].T.values
    beta, *_ = np.linalg.lstsq(Xt, Yt, rcond=None)
    resid = Yt - Xt @ beta; dof = len(train) - Xt.shape[1]
    se = np.sqrt((resid ** 2).sum(axis=0) / dof * np.linalg.inv(Xt.T @ Xt)[1, 1])
    tval = beta[1] / se; p = 2 * stats.t.sf(np.abs(tval), dof)
    return [g for g, b, pp in zip(hyp_in, beta[1], p) if b >= 1 and pp < 0.001]
def score_one(sample, train, genes):
    base = Lh.loc[genes, [s for s in train if metadata.loc[s, "pe_status"] == "no_PE"]].mean(axis=1)
    return (Lh.loc[genes, sample] - base).mean()
cv_score, n_sel = {}, {}
for s in preterm:
    train = [x for x in preterm if x != s]; g = select_genes(train)
    cv_score[s] = score_one(s, train, g); n_sel[s] = len(g)
genes_all = select_genes(list(preterm))
for s in metadata.index[metadata.group == "TB"]:
    cv_score[s] = score_one(s, list(preterm), genes_all); n_sel[s] = len(genes_all)
metadata["cv_hypoxia_score"] = pd.Series(cv_score)
print(f"Genes chosen per fold: median {np.median(list(n_sel.values())):.0f} (range {min(n_sel.values())} to {max(n_sel.values())}). "
      f"From all 25: {', '.join(genes_all)}")
fit = smf.ols("cv_hypoxia_score ~ C(pe_status, Treatment('no_PE')) + C(fgr_status) + C(fetal_sex) + ga_delivery_weeks",
              metadata.loc[preterm]).fit()
t = [x for x in fit.params.index if x.startswith("C(pe_status")][0]; ci = fit.conf_int().loc[t]
headline_stats = dict(effect=fit.params[t], lo=ci[0], hi=ci[1], p=fit.pvalues[t], slope=fit.params["ga_delivery_weeks"],
                      slope_p=fit.pvalues["ga_delivery_weeks"])
print("Cross-validated score, PE effect: %.2f log2 (%.1fx, 95%% CI %.1fx to %.1fx), p = %.1e; delivery-week slope %.3f per week (p = %.2f)" %
      (headline_stats["effect"], 2 ** headline_stats["effect"], 2 ** headline_stats["lo"], 2 ** headline_stats["hi"],
       headline_stats["p"], headline_stats["slope"], headline_stats["slope_p"]))
nonpe_pre = preterm[metadata.loc[preterm, "pe_status"] == "no_PE"]
print("Term births vs preterm non-PE: Mann-Whitney p = %.2f" %
      stats.mannwhitneyu(metadata.loc[metadata.group == "TB", "cv_hypoxia_score"], metadata.loc[nonpe_pre, "cv_hypoxia_score"]).pvalue)
metadata[["group", "ga_delivery_weeks", "cv_hypoxia_score"]].to_csv(RESULTS_DIR / "per_sample_cv_score.csv")

# %%
fig, ax = plt.subplots(figsize=(9, 5.8))
sc = metadata["cv_hypoxia_score"]
for g, st in GROUP_STYLE.items():
    idx = metadata.index[metadata.group == g]
    ax.scatter(metadata.loc[idx, "ga_delivery_weeks"], 2 ** sc[idx], s=95, facecolor=st["color"], marker=st["marker"],
               edgecolor="white", linewidth=1.3, label=f'{st["label"]} ({len(idx)})', zorder=3)
ax.axhline(1, color=MUTED, lw=0.9, ls="--", zorder=1)
ax.axvspan(34, 38.8, color="#f2f1ee", zorder=0)
ax.text(36.4, 6, "no deliveries\nin this range", ha="center", fontsize=8.5, color=MUTED)
ax.text(39.6, 1.07, "average of preterm\nbirths without\npreeclampsia", fontsize=8.5, color=MUTED, ha="right", va="bottom")
ax.set_yscale("log", base=2)
ax.set_yticks([0.5, 1, 2, 4, 8]); ax.set_yticklabels(["½×", "1×", "2×", "4×", "8×"])
ax.set_xlim(25, 40)
ax.set_xlabel("Weeks of pregnancy at delivery")
ax.set_ylabel("Activity of low-oxygen response genes\n(relative to preterm births without preeclampsia)")
ax.set_title("Hypoxia-Response Score by PE Status and Delivery Week", loc="left", fontsize=14)
ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, scatterpoints=1)
fig.text(0.01, -0.03, "Each point is one placenta (n = 30). Genes for each placenta's score were chosen without that placenta "
         "(cross-validated). Data: GSE203507.", fontsize=8, color=MUTED, ha="left")
fig.tight_layout(); save(fig, "headline.png")

# %%
fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True)
axes = axes.ravel()
for ax, gene in zip(axes, ["FLT1", "ENG", "LEP", "HLA-G"]):
    for g, st in GROUP_STYLE.items():
        idx = metadata.index[metadata.group == g]
        ax.scatter(metadata.loc[idx, "ga_delivery_weeks"], log_cpm_all.loc[gene, idx], s=55, facecolor=st["color"],
                   marker=st["marker"], edgecolor="white", linewidth=1, label=st["label"])
    if gene == "HLA-G":
        ax.set_title("HLA-G by Group and Delivery Week")
    else:
        ax.set_title(f"{gene}: {key.loc[gene, 'fold']:.1f}× Higher in Preeclampsia")
    ax.set_xlabel("Weeks of pregnancy at delivery")
    ax.set_ylabel("Gene activity (log2 CPM)")
print("HLA-G: same-week PE effect log2FC = %.2f (padj %.2f); term vs preterm controls log2FC = %.2f (padj %.3f)" %
      (primary.loc["HLA-G", "log2FoldChange"], primary.loc["HLA-G", "padj"], ga_sig.loc["HLA-G", "log2FoldChange"], ga_sig.loc["HLA-G", "padj"]))
fig.tight_layout(rect=(0, 0.06, 1, 1))
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=9, loc="lower center", ncol=3, scatterpoints=1)
save(fig, "fig12_key_genes.png")

# %%
# Do the clinical biomarker genes (not in the Hallmark set) move with the hypoxia program? Within PE only.
pe_idx = metadata.index[metadata.pe_status == "PE"]
for gene in ["FLT1", "ENG", "LEP", "INHBA", "FSTL3", "HTRA4"]:
    r = stats.spearmanr(log_cpm_all.loc[gene, pe_idx], metadata.loc[pe_idx, "cv_hypoxia_score"])
    print(f"{gene}: rho with the hypoxia score within the 15 PE placentas = {r[0]:.2f} (p = {r[1]:.3f})")
for v in ["sbp_max_mmhg", "dbp_max_mmhg", "ga_delivery_weeks"]:
    r = stats.spearmanr(metadata.loc[pe_idx, "cv_hypoxia_score"], metadata.loc[pe_idx, v])
    print(f"Within PE (n=15): score vs {v}: rho = {r[0]:.2f}, p = {r[1]:.2f}")
display(metadata.loc[nonpe_pre, ["group", "ga_delivery_weeks", "sbp_max_mmhg", "dbp_max_mmhg", "qc_flags", "cv_hypoxia_score"]]
        .sort_values("cv_hypoxia_score", ascending=False).round(2))

# %% [markdown]
# ## 12. Summary of Answers
#
# Printed from the objects above so the numbers stay in sync with the article.

# %%
print("Q6. How much of the PE vs control difference is delivery timing?")
print(f"  - PC1 ({var[0]:.0%} of variation) tracks PE (p = {assoc.loc['Preeclampsia','PC1']:.0e}), not delivery week (p = {assoc.loc['Weeks at delivery','PC1']:.2f}).")
print(f"  - Median share of each gene's PE difference explained by timing: {q6['term indicator']['median']:.0%} (term indicator), "
      f"{q6['weeks']['median']:.0%} (weeks); {q6['term indicator']['over_half']} of {len(share)} genes lose more than half.")
print(f"  - Partial R2 in PE genes: PE {vp_summary.loc['Genes passing (no-timing model)','Preeclampsia']:.0%} vs weeks "
      f"{vp_summary.loc['Genes passing (no-timing model)','Weeks at delivery']:.0%} (medians).")
print(f"  - Timing-matched comparison: {len(hits)} genes; overlap with the delivery-timing signature: {len(set(h_ga.index) & set(hits.index))}.")
print(f"  - Timing signature {len(h_ga)} genes vs random control splits {min(noise_draws)} to {max(noise_draws)} and 5 PE vs 5 preterm "
      f"{min(disease_draws)} to {max(disease_draws)}.")
print(f"  - Worst case (term controls only): median timing part {worst_case_median:.0%} of each PE-vs-term difference.")
print("Q9. Broad or concentrated?")
print(f"  - Median tau: higher-in-PE {spec.iloc[0].median_tau:.2f} vs matched random {spec.iloc[0].null_median_tau:.2f} (p = {spec.iloc[0].empirical_p:.3f}); "
      f"all expressed {spec.iloc[3].median_tau:.2f}.")
print(f"  - {n_up_t} of {n_up} higher-in-PE genes peak in trophoblast in healthy placenta vs {bg_t:.0%} background "
      f"and {timing_t} of {timing_n} delivery-timing genes.")
stb_late_row = conc.query("onset=='late' and cell_class=='syncytiotrophoblast'").iloc[0]
print(f"  - Within syncytiotrophoblast, late onset (timing-matched): {stb_late_row.frac_higher_in_pe:.0%} higher in PE, median +{stb_late_row.median_log2_diff:.2f} log2 (null p = {stb_late_row.null_p:.3f}).")
print("Headline: cross-validated hypoxia score %.1fx higher in PE (95%% CI %.1f to %.1f), slope with delivery week p = %.2f." %
      (2 ** headline_stats["effect"], 2 ** headline_stats["lo"], 2 ** headline_stats["hi"], headline_stats["slope_p"]))

# %%
import pydeseq2, sklearn, scipy, statsmodels
print("python", sys.version.split()[0], "| numpy", np.__version__, "| pandas", pd.__version__, "| pydeseq2", pydeseq2.__version__,
      "| scikit-learn", sklearn.__version__, "| scipy", scipy.__version__, "| statsmodels", statsmodels.__version__)
