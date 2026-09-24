#!/usr/bin/env python3
"""Plot tested human miRNAs using the source paper's display thresholds.

The candidate panel is selected from the corrected main DESeq2 results using
unadjusted p < 0.05 and at least a twofold change (absolute log2FC >= 1).
This is exploratory; FDR is assessed over all finite tested features.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error(f"Output exists: {a.out}")
    a.out.mkdir(parents=True)
    de = pd.read_csv(a.models / "main_ge2000/deseq2_results.csv")
    ed = pd.read_csv(a.models / "main_ge2000/edger_TMMwsp_results.csv")
    full = pd.read_csv(a.models / "model_summary.csv").set_index("strategy")
    finite = de[np.isfinite(de.pvalue) & np.isfinite(de.log2FoldChange)].copy()
    selected = (finite.pvalue < .05) & (finite.log2FoldChange.abs() >= 1)
    nominal = finite[selected].sort_values("pvalue").copy()
    assert len(de) == 244 and len(finite) == 230 and len(nominal) == 15
    assert not (nominal.padj < .05).any() and not (finite.padj < .05).any()
    assert full.loc["main_ge2000", "deseq2_FDR05"] == 0
    assert full.loc["main_ge2000", "edger_FDR05"] == 0
    nominal = nominal.merge(ed[["miRNA", "logFC", "PValue", "FDR"]],
                            on="miRNA", how="left", validate="one_to_one")
    nominal["direction"] = np.where(nominal.log2FoldChange > 0,
                                     "higher in cases", "lower in cases")
    nominal["all22_log2FoldChange"] = nominal.miRNA.map(
        pd.read_csv(a.models / "all22/deseq2_results.csv").set_index("miRNA").log2FoldChange)
    nominal[["miRNA", "direction", "baseMean", "log2FoldChange", "wald_lower95",
             "wald_upper95", "pvalue", "padj", "logFC", "PValue", "FDR",
             "all22_log2FoldChange"]].rename(columns={
                 "pvalue": "DESeq2_p", "padj": "DESeq2_FDR",
                 "logFC": "edgeR_log2FC", "PValue": "edgeR_p", "FDR": "edgeR_FDR"
             }).to_csv(a.out / "Table_2_revised_nominal_human_candidates.csv", index=False)

    up, down, other = "#b04a36", "#31835a", "#b8bec4"
    colors = np.where(selected,
                      np.where(finite.log2FoldChange > 0, up, down), other)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.3, 7.2),
                                 gridspec_kw={"width_ratios": [1.05, 1.2]},
                                 layout="constrained")
    ax.scatter(finite.log2FoldChange, -np.log10(finite.pvalue), c=colors,
               s=24, alpha=.86, linewidths=0)
    ax.axhline(-np.log10(.05), color="#666666", lw=.9, ls="--")
    ax.axvline(0, color="#777777", lw=.7)
    for threshold in (-1, 1):
        ax.axvline(threshold, color="#999999", lw=.8, ls="--")
    ax.set(xlabel="Case versus control log2 fold change",
           ylabel="−log10(unadjusted p)",
           title="A  All 230 miRNAs with finite tests")
    ax.text(.03, .97, "15 with p < 0.05 and ≥2-fold change\n0 with q < 0.05",
            transform=ax.transAxes, va="top", ha="left",
            bbox={"facecolor": "white", "edgecolor": "#dddddd", "pad": 6})
    ax.scatter([], [], c=up, s=32, label="Higher in cases; display threshold met")
    ax.scatter([], [], c=down, s=32, label="Lower in cases; display threshold met")
    ax.scatter([], [], c=other, s=32, label="Other finite tests")
    ax.legend(loc="lower left", frameon=False, fontsize=8)

    y = np.arange(len(nominal))
    for i, row in enumerate(nominal.itertuples()):
        color = up if row.log2FoldChange > 0 else down
        bx.errorbar(row.log2FoldChange, i,
                    xerr=np.array([[row.log2FoldChange - row.wald_lower95],
                                   [row.wald_upper95 - row.log2FoldChange]]),
                    fmt="o", color=color, ecolor=color, capsize=2, markersize=5)
    bx.axvline(0, color="#777777", lw=.8)
    bx.set(yticks=y, yticklabels=nominal.miRNA,
           xlabel="Case versus control log2 fold change (95% Wald interval)",
           title="B  Revised candidates meeting both display thresholds")
    bx.invert_yaxis()
    bx.grid(axis="x", color="#e2e2e2", lw=.5)
    lo = min(-4.8, nominal.wald_lower95.min() - .4)
    hi = max(5.5, nominal.wald_upper95.max() + .3)
    bx.set_xlim(lo, hi + 2.2)
    for i, row in enumerate(nominal.itertuples()):
        bx.text(hi + .22, i, f"q={row.padj:.3f}", va="center", fontsize=8,
                color="#555555")
    fig.suptitle("Human miRNAs in corrected-count analysis",
                 fontsize=12, fontweight="bold")
    fig.savefig(a.out / "Figure_1_human_miRNA.png", dpi=220)
    fig.savefig(a.out / "Figure_1_human_miRNA.svg")
    plt.close(fig)
    summary = {"modeled_features": len(de), "finite_DESeq2_p": len(finite),
               "original_style_p_below_0_05_and_abs_log2FC_ge_1": len(nominal),
               "nominal_higher_cases": int((nominal.log2FoldChange > 0).sum()),
               "nominal_lower_cases": int((nominal.log2FoldChange < 0).sum()),
               "FDR_below_0_05": int((finite.padj < .05).sum()),
               "minimum_FDR": float(finite.padj.min()),
               "selection_note": "Source-paper-style display rule (unadjusted p<0.05 and >=2-fold change) applied to the revised DESeq2 model; not a replication of its CLC workflow or a validated differential signature."}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
