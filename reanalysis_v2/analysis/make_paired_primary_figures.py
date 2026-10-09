"""Plot both primary modalities in the same 16 recipients."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sankey_figure import draw_sankey


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paired_primary_2026-10-09/figures"
OUT.mkdir(parents=True, exist_ok=True)
MODEL = ROOT / "joint_mirna_models_2026-10-02_v2/main_ge2000"

de = pd.read_csv(MODEL / "deseq2_results.csv")
host = de[de.feature_type.eq("human")]
finite = host[np.isfinite(host.pvalue) & np.isfinite(host.log2FoldChange)].copy()
candidates = finite[(finite.pvalue < .05) & (finite.log2FoldChange.abs() >= 1)].sort_values("pvalue")
assert len(finite) == 230 and len(candidates) == 15
viral = de[de.feature_type.eq("viral")].set_index("miRNA")
assert np.isfinite(viral.loc["bkv-miR-B1-3p", "pvalue"])
assert np.isnan(viral.loc["bkv-miR-B1-5p", "pvalue"])

metadata = pd.read_csv(ROOT / "verified_data/metadata/bkv_sample_metadata.csv").set_index("sample")
classes = pd.read_csv(ROOT / "results/tables/class_profiles_all22.csv", index_col=0)
primary_samples = pd.read_csv(MODEL / "fixed_host_normalization.csv")["sample"].tolist()
controls = [s for s in primary_samples if metadata.loc[s, "group"] == "control"]
cases = [s for s in primary_samples if metadata.loc[s, "group"] == "treatment"]
assert len(controls) == 7 and len(cases) == 9
order = sorted(controls, key=lambda s: int(s[1:])) + sorted(cases, key=lambda s: int(s[1:]))

plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "svg.fonttype": "none",
})
up, down, other, virus = "#ac4a37", "#267c54", "#bac1c7", "#713d9c"
fig, (ax, bx) = plt.subplots(
    1, 2, figsize=(18, 7.5), gridspec_kw={"width_ratios": [1, 1.85]}
)
fig.subplots_adjust(left=.065, right=.985, top=.90, bottom=.19, wspace=.13)

# A: the exact primary joint DESeq2 model used in the manuscript.
chosen = (finite.pvalue < .05) & (finite.log2FoldChange.abs() >= 1)
colors = np.where(chosen, np.where(finite.log2FoldChange > 0, up, down), other)
ax.scatter(finite.log2FoldChange, -np.log10(finite.pvalue), c=colors,
           s=26, alpha=.86, linewidths=0, zorder=2)
v3, v5 = viral.loc["bkv-miR-B1-3p"], viral.loc["bkv-miR-B1-5p"]
ax.scatter([v3.log2FoldChange], [-np.log10(v3.pvalue)], marker="*", s=160,
           c=virus, edgecolors="white", linewidths=.6, zorder=5)
ax.axhspan(-.72, -.04, color="#f4f0f7", zorder=0)
ax.scatter([v5.log2FoldChange], [-.39], marker="D", s=85,
           facecolors="none", edgecolors=virus, linewidths=1.6, zorder=5)
ax.text(v5.log2FoldChange-.12, -.39, "5p: p unavailable", ha="right",
        va="center", color=virus, fontsize=8)
ax.axhline(-np.log10(.05), color="#777", lw=.8, ls="--")
ax.axhline(0, color="#999", lw=.6)
ax.axvline(0, color="#777", lw=.7)
for boundary in (-1, 1):
    ax.axvline(boundary, color="#aaa", lw=.7, ls="--")
ax.set_ylim(-.76, 5)
ax.set_xlim(-4.1, 5.3)
ax.set_yticks(range(6))
ax.set(xlabel="Case versus control log2 fold change",
       ylabel="−log10(unadjusted p)",
       title="A  Human and viral miRNAs (9 cases, 7 controls)")
ax.text(.98, .98, "15 human display candidates; 0 q < 0.05\n"
        "Shared viral 3p q = 0.057", transform=ax.transAxes, ha="right",
        va="top", fontsize=8,
        bbox={"facecolor":"white", "edgecolor":"#ddd", "pad":5})
for color, label, marker in [
    (up, "Human higher in cases; nominal rule met", "o"),
    (down, "Human lower in cases; nominal rule met", "o"),
    (other, "Other tested human miRNAs", "o"),
    (virus, "BKPyV/JCPyV-shared 3p", "*"),
]:
    ax.scatter([], [], c=color, s=48 if marker == "*" else 27,
               marker=marker, label=label)
ax.scatter([], [], marker="D", s=35, facecolors="none", edgecolors=virus,
           label="BKPyV 5p: Cook's-filtered")
ax.legend(loc="upper left", bbox_to_anchor=(.0, .99), frameon=False, fontsize=7)

# B: the same 16 individual bacterial profiles, then equal-weight group means.
# All 22 profiles remain available for the secondary sensitivity analysis.
major = classes[order].mean(axis=1).nlargest(7).index.tolist()
display = classes.loc[major, order].copy()
display.loc["Other classes"] = classes.drop(major).loc[:, order].sum()
assert np.allclose(display.sum(axis=0), 100, atol=.02)
means = pd.DataFrame({
    "Control mean": display[sorted(controls, key=lambda s: int(s[1:]))].mean(axis=1),
    "Case mean": display[sorted(cases, key=lambda s: int(s[1:]))].mean(axis=1),
})
plotted = pd.concat([display, means], axis=1)
x = np.r_[np.arange(7), np.arange(8,17), [18,19]]
palette = plt.get_cmap("tab10")(np.arange(len(plotted)))
bottom = np.zeros(len(x))
for (_, values), color in zip(plotted.iterrows(), palette):
    bx.bar(x, values.to_numpy(), bottom=bottom, width=.81, color=color)
    bottom += values.to_numpy()
bx.axvline(7, color="#555", lw=.7)
bx.axvline(17.2, color="#555", lw=.7)
bx.set_ylim(0, 100)
bx.set_xlim(-.8, 19.8)
bx.set_xticks(x)
bx.set_xticklabels([*order, "Ctrl", "Case"], rotation=90, fontsize=7)
bx.set(ylabel="Bacterial relative abundance (%)",
       title="B  MetaPhlAn bacterial classes (same 16 recipients)")
bx.text(3, -12, "Blood negative", ha="center", va="top", fontsize=8)
bx.text(12, -12, "Blood DNAemia", ha="center", va="top", fontsize=8)
bx.text(18.5, -12, "Means", ha="center", va="top", fontsize=8)
handles = [plt.Rectangle((0,0),1,1, color=color) for color in palette]
fig.legend(handles, list(plotted.index), loc="lower center", ncol=4,
           frameon=False, bbox_to_anchor=(.69, .01), fontsize=8)
fig.suptitle("Urinary miRNAs and bacterial composition by blood BKPyV DNAemia",
             fontsize=13, fontweight="bold")
for ext in ("png", "svg"):
    fig.savefig(OUT / f"Figure_1_miRNA_MetaPhlAn.{ext}",
                dpi=240 if ext == "png" else None)
plt.close(fig)

# Figure 2 uses the same 16 recipients and 15 human features as Figure 1A.
sf = pd.read_csv(MODEL / "fixed_host_normalization.csv").set_index("sample").host_size_factor
paired = sf.index.tolist()
paired_controls = [s for s in paired if metadata.loc[s, "group"] == "control"]
paired_cases = [s for s in paired if metadata.loc[s, "group"] == "treatment"]
assert len(paired_controls) == 7 and len(paired_cases) == 9
human = pd.read_csv(ROOT / "inputs/recount/human_unique_counts.csv", index_col=0)
virus_counts = pd.read_csv(ROOT / "inputs/recount/viral_historical_counts.csv", index_col=0)
features = candidates.miRNA.tolist() + ["bkv-miR-B1-3p", "bkv-miR-B1-5p"]
normalized = np.log2(pd.concat([human, virus_counts])[paired].div(sf, axis=1) + 1).loc[features]
draw_sankey(classes[paired], normalized, paired_controls, paired_cases,
            OUT, OUT, figure_number=2, right_label="Current miRNA display candidates")
assert "hsa-miR-320b" not in features

(OUT / "Figure_1_caption.txt").write_text(
    "Figure 1. Urinary miRNAs and MetaPhlAn bacterial class profiles by blood BKPyV "
    "DNAemia. A, Volcano plot from the primary joint DESeq2 model of nine cases "
    "and seven controls. Human points with unadjusted P<0.05 and |log2 fold change| "
    "≥1 are colored by direction; no human or viral feature had q<0.05. The star "
    "is the sequence shared by BKPyV B1-3p and JCPyV J1-3p (q=0.057). The hollow "
    "diamond shows the BKPyV 5p estimated effect in an untested band because "
    "Cook's-distance filtering removed its primary P value. B, MetaPhlAn "
    "class-level relative abundances for the same 16 recipients and means of the "
    "individual control and case profiles. Other classes comprise all classes "
    "outside the seven most abundant in these 16. Each stacked bar totals 100%; "
    "the mean bars were calculated across recipients, not pooled reads. "
    "These profiles describe composition, not bacterial load or viable organisms.\n"
)
(OUT / "Figure_2_caption.txt").write_text(
    "Figure 2. Descriptive Sankey connecting bacterial class profiles and miRNA "
    "patterns by sample identity in the paired nine-case/seven-control subset. "
    "Left ribbons encode bacterial class relative percentages; classes below "
    "1% in every displayed sample are grouped as Other, and dashed guides show "
    "tiny connections without encoding abundance. Right ribbons encode absolute "
    "deviations of log2(count/human DESeq2 size factor + 1) from the control mean "
    "for the 15 Figure 1A human display candidates and two viral sequences; "
    "deviations <0.5 are omitted. The two sides have independent units and "
    "width scales. Human above/below-control deviations are red/green, and viral "
    "above/below deviations are dark/light blue. The asterisk marks the "
    "BKPyV/JCPyV-shared 3p sequence. Links identify the same sample on both "
    "sides and do not represent correlations, interactions, or causal effects.\n"
)
print(OUT)
