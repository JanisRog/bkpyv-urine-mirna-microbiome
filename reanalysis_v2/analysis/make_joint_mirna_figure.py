#!/usr/bin/env python3
"""Figure 1: human and viral miRNAs from the same joint DESeq2 model."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"Output directory exists: {args.out}")
    args.out.mkdir(parents=True)
    de = pd.read_csv(args.models / "main_ge2000/deseq2_results.csv")
    assert len(de) == 246 and (de.feature_type == "human").sum() == 244
    host = de[de.feature_type.eq("human")].copy()
    viral = de[de.feature_type.eq("viral")].set_index("miRNA")
    finite_host = host[np.isfinite(host.pvalue) & np.isfinite(host.log2FoldChange)].copy()
    candidates = finite_host[(finite_host.pvalue < 0.05) &
                             (finite_host.log2FoldChange.abs() >= 1)].sort_values("pvalue")
    assert len(finite_host) == 230 and len(candidates) == 15
    assert np.isfinite(viral.loc["bkv-miR-B1-3p", "pvalue"])
    assert np.isnan(viral.loc["bkv-miR-B1-5p", "pvalue"])
    assert not (de.padj < .05).any()
    up, down, other, virus = "#ac4a37", "#267c54", "#bac1c7", "#713d9c"
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15.6, 8.2),
                                 gridspec_kw={"width_ratios": [1, 1.2]},
                                 layout="constrained")
    selected = (finite_host.pvalue < .05) & (finite_host.log2FoldChange.abs() >= 1)
    colors = np.where(selected,
                      np.where(finite_host.log2FoldChange > 0, up, down), other)
    ax.scatter(finite_host.log2FoldChange, -np.log10(finite_host.pvalue),
               c=colors, s=25, alpha=.86, linewidths=0, zorder=2)
    v3 = viral.loc["bkv-miR-B1-3p"]
    v5 = viral.loc["bkv-miR-B1-5p"]
    ax.scatter([v3.log2FoldChange],[-np.log10(v3.pvalue)], s=160, marker="*",
               c=virus, edgecolors="white",linewidths=.7,zorder=5,
               label="BKPyV/JCPyV shared 3p")
    # A missing p-value cannot be assigned a volcano ordinate. Show its
    # estimated effect in a separate, explicitly untested band.
    ax.axhspan(-.72,-.04,color="#f4f0f7",zorder=0)
    ax.scatter([v5.log2FoldChange],[-.39],s=86,marker="D",facecolors="none",
               edgecolors=virus,linewidths=1.7,zorder=5,
               label="BKPyV 5p: Cook's-filtered")
    ax.text(v5.log2FoldChange-.13,-.38,"5p: p unavailable",ha="right",va="center",
            color=virus,fontsize=8)
    ax.axhline(-np.log10(.05),color="#777",lw=.85,ls="--")
    ax.axhline(0,color="#999",lw=.6)
    ax.axvline(0,color="#777",lw=.7)
    for x in (-1,1):
        ax.axvline(x,color="#aaa",lw=.7,ls="--")
    ax.set(xlabel="Case versus control log2 fold change",
           ylabel="−log10(unadjusted p)",
           title="A  Joint human and viral count model")
    ax.set_ylim(bottom=-.76)
    ax.set_yticks([0,1,2,3,4,5])
    ax.text(.98,.98,"15 human display candidates; 0 q < 0.05\n"
            "Viral 3p q = 0.057; 5p test filtered",
            ha="right",va="top",transform=ax.transAxes,fontsize=8,
            bbox={"facecolor":"white","edgecolor":"#ddd","pad":6})
    ax.scatter([],[],c=up,s=28,label="Human: higher in cases, display rule met")
    ax.scatter([],[],c=down,s=28,label="Human: lower in cases, display rule met")
    ax.scatter([],[],c=other,s=28,label="Other human finite tests")
    ax.legend(loc="lower left",bbox_to_anchor=(.0,.06),frameon=False,fontsize=7.5)

    rows = pd.concat([candidates.set_index("miRNA"), viral],axis=0)
    y_positions = list(range(len(candidates))) + [len(candidates)+.8,
                                              len(candidates)+1.8]
    labels = list(candidates.miRNA) + ["BKPyV/JCPyV shared 3p", "BKPyV 5p"]
    for i, (name,row) in enumerate(rows.iterrows()):
        yy = y_positions[i]
        if name == "bkv-miR-B1-5p":
            bx.scatter([row.log2FoldChange],[yy],marker="D",s=55,
                       facecolors="none",edgecolors=virus,linewidths=1.7,zorder=4)
        else:
            color = virus if row.feature_type == "viral" else (up if row.log2FoldChange > 0 else down)
            bx.errorbar(row.log2FoldChange,yy,
                        xerr=np.array([[row.log2FoldChange-row.wald_lower95],
                                       [row.wald_upper95-row.log2FoldChange]]),
                        fmt="*" if row.feature_type=="viral" else "o",
                        color=color,ecolor=color,capsize=2,markersize=9 if row.feature_type=="viral" else 5)
    bx.axvline(0,color="#777",lw=.8)
    bx.axhline(len(candidates)-.1,color="#d9d3de",lw=.8)
    bx.set(yticks=y_positions,yticklabels=labels,
           xlabel="Case versus control log2 fold change (95% Wald interval where shown)",
           title="B  Human display candidates and prespecified viral sequences")
    bx.invert_yaxis()
    bx.grid(axis="x",color="#e4e4e4",lw=.5)
    lo = min(-4.8, np.nanmin(rows.wald_lower95)-.4)
    hi = max(6.0, np.nanmax(rows.wald_upper95)+.3)
    bx.set_xlim(lo,hi+2.8)
    for i, (_,row) in enumerate(rows.iterrows()):
        label = "p unavailable" if pd.isna(row.pvalue) else f"q={row.padj:.3f}"
        bx.text(hi+.22,y_positions[i],label,va="center",fontsize=7.5,
                color=virus if row.feature_type=="viral" else "#555")
    fig.suptitle("Human and viral urinary miRNAs in the 9-case/7-control analysis",
                 fontsize=12,fontweight="bold")
    fig.savefig(args.out / "Figure_1_joint_human_viral_miRNA.png",dpi=240)
    fig.savefig(args.out / "Figure_1_joint_human_viral_miRNA.svg")
    plt.close(fig)

    plot_data = rows.reset_index(names="miRNA").copy()
    plot_data.loc[plot_data.miRNA.eq("bkv-miR-B1-5p"),
                  ["wald_lower95","wald_upper95"]] = np.nan
    plot_data[
        ["miRNA","feature_type","baseMean","log2FoldChange","wald_lower95",
         "wald_upper95","pvalue","padj","max_Cooks"]
    ].to_csv(args.out / "Figure_1_plot_data.csv",index=False)
    summary = {
        "model":"Joint DESeq2 negative-binomial model with fixed human-derived size factors",
        "samples":"9 cases and 7 controls", "features":246,
        "human_finite_tests":230,"viral_finite_tests":1,
        "human_display_candidates":15,"q_below_0_05":0,
        "viral_3p_q":float(v3.padj),
        "viral_5p_note":"Cook's-distance filtered in primary model; effect shown without an inferential interval or p-value."
    }
    (args.out / "summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (args.out / "Figure_1_caption.txt").write_text(
        "Figure 1. Human and viral urinary miRNAs analyzed in a common DESeq2 "
        "negative-binomial model for nine blood-BKPyV-DNAemia cases and seven controls. "
        "Size factors were estimated from human miRNAs before the viral rows were "
        "added. A, All 230 human features with finite tests and the BKPyV/JCPyV-shared "
        "3p sequence are placed by unadjusted Wald p-value; colors identify human "
        "features meeting the prespecified descriptive display rule (p<0.05 and "
        "|log2 fold change|>=1). The BKPyV 5p estimate is placed in a separate "
        "untested band because Cook's-distance filtering removed its p-value. "
        "B, The 15 human display candidates and both prespecified viral sequences; "
        "horizontal intervals are unshrunk 95% Wald intervals for features with "
        "valid primary tests. BKPyV 5p is shown as a hollow diamond without an "
        "inferential interval. q values are Benjamini-Hochberg adjusted across "
        "231 finite tests; none is below 0.05. The shared 3p sequence cannot be "
        "assigned specifically to BKPyV rather than JCPyV from this assay.\n"
    )
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
