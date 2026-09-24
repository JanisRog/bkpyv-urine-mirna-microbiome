#!/usr/bin/env python3
"""Make sample-level DNA/RNA evidence tables and readable final revision figures.

This is an exploratory integration added after the prespecified table analyses.
Counts from different libraries are deliberately displayed on separate axes.
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig, directory, stem):
    fig.savefig(directory / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(directory / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--alignment", type=Path, required=True)
    p.add_argument("--kraken", type=Path, required=True)
    p.add_argument("--recount", type=Path, required=True)
    p.add_argument("--qc", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error(f"Output already exists: {a.out}")
    a.out.mkdir(parents=True)
    (a.out / "figures").mkdir()
    (a.out / "tables").mkdir()
    figdir, tabdir = a.out / "figures", a.out / "tables"
    meta = pd.read_csv(a.root / "verified_data/metadata/bkv_sample_metadata.csv")
    meta["group"] = meta.group.replace({"treatment": "case"})
    assert len(meta) == 22 and meta["sample"].is_unique and set(meta.group) == {"case", "control"}
    align = pd.read_csv(a.alignment, sep="\t")
    align["group"] = align.group.replace({"treatment": "case"})
    assert len(align) == 44 and not align.duplicated(["sample", "reference"]).any()
    assert set(align.reference) == {"NC_001538.1", "NC_001699.1"}
    assert set(align["sample"]) == set(meta["sample"])
    kraken = pd.read_csv(a.kraken, sep="\t")
    assert len(kraken) == 22 and kraken.accession.is_unique
    assert set(kraken.accession) == set(align.accession)
    viral = pd.read_csv(a.recount / "viral_historical_counts.csv", index_col=0)
    assert {"bkv-miR-B1-3p", "bkv-miR-B1-5p"} <= set(viral.index)
    qc = pd.read_csv(a.qc)
    assert set(qc["sample"]) == set(meta["sample"])
    order = (meta.query("group == 'control'")["sample"].tolist() +
             meta.query("group == 'case'")["sample"].tolist())
    al = align.set_index(["sample", "reference"])
    ka = kraken.set_index("sample")
    qua = qc.set_index("sample")
    ma = meta.set_index("sample")
    records = []
    for sample in order:
        bk, jc = al.loc[(sample, "NC_001538.1")], al.loc[(sample, "NC_001699.1")]
        k = ka.loc[sample]
        assert bk.accession == jc.accession == k.accession
        assert bk.group == jc.group == ma.loc[sample, "group"]
        records.append({
            "sample": sample, "accession": bk.accession, "group": bk.group,
            "host_filtered_dna_pairs": int(k.total_pairs),
            "kraken_bkpyv_pairs": int(k.bkpyv_pairs),
            "kraken_jcpyv_pairs": int(k.jcpyv_pairs),
            "bkpyv_mapq30_first_mates": int(bk.mapq30_first_mates),
            "jcpyv_mapq30_first_mates": int(jc.mapq30_first_mates),
            "bkpyv_mapq30_first_mates_pct_host_filtered_pairs": 100 * int(bk.mapq30_first_mates) / int(k.total_pairs),
            "jcpyv_mapq30_first_mates_pct_host_filtered_pairs": 100 * int(jc.mapq30_first_mates) / int(k.total_pairs),
            "bkpyv_breadth_1x_pct": bk.breadth_1x_pct_mapq30_baseq20,
            "bkpyv_breadth_10x_pct": bk.breadth_10x_pct_mapq30_baseq20,
            "bkpyv_breadth_100x_pct": bk.breadth_100x_pct_mapq30_baseq20,
            "jcpyv_breadth_1x_pct": jc.breadth_1x_pct_mapq30_baseq20,
            "jcpyv_breadth_10x_pct": jc.breadth_10x_pct_mapq30_baseq20,
            "jcpyv_breadth_100x_pct": jc.breadth_100x_pct_mapq30_baseq20,
            "bkpyv_median_depth": bk.median_depth_mapq30_baseq20,
            "jcpyv_median_depth": jc.median_depth_mapq30_baseq20,
            "shared_bk_jc_3p_mirna_raw_reads": int(viral.loc["bkv-miR-B1-3p", sample]),
            "bkpyv_5p_mirna_raw_reads": int(viral.loc["bkv-miR-B1-5p", sample]),
            "smallrna_retained_reads": int(qua.loc[sample, "retained_reads"]),
            "kraken_fungi_pairs": int(k.fungi_pairs),
            "kraken_human_pairs": int(k.human_pairs),
            "kraken_fungi_pct_classifier_input": 100 * k.fungi_pairs / k.total_pairs,
            "kraken_human_pct_classifier_input": 100 * k.human_pairs / k.total_pairs,
        })
    df = pd.DataFrame(records)
    assert (df.bkpyv_breadth_100x_pct > 95).all()
    df.to_csv(tabdir / "sample_level_viral_and_qc.csv", index=False)
    checksums = {str(x): sha256(x) for x in [a.alignment, a.kraken,
                 a.recount / "viral_historical_counts.csv", a.qc,
                 a.root / "verified_data/metadata/bkv_sample_metadata.csv"]}
    (a.out / "input_sha256.json").write_text(json.dumps(checksums, indent=2) + "\n")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False})
    colors = np.where(df.group.eq("case"), "#a94a34", "#236c91")
    x = np.arange(len(df))
    fig, axs = plt.subplots(4, 1, figsize=(12, 9), sharex=True,
                            layout="constrained", gridspec_kw={"height_ratios": [1, 1, 1, 1.1]})
    for ax in axs:
        ax.axvspan(10.5, 21.5, color="#a94a34", alpha=.035)
        ax.axvline(10.5, color="#777777", linewidth=.8)
        ax.grid(axis="y", color="#dddddd", linewidth=.5)
    axs[0].scatter(x, df.bkpyv_mapq30_first_mates + 1, c=colors, s=38)
    axs[0].set(yscale="log", ylabel="BKPyV high-MAPQ\nfirst mates + 1")
    axs[1].scatter(x, df.jcpyv_mapq30_first_mates + 1, c=colors, s=38)
    axs[1].set(yscale="log", ylabel="JCPyV high-MAPQ\nfirst mates + 1")
    axs[1].annotate("S14: deep JC coverage", (df.index[df["sample"].eq("S14")][0],
                    float(df.loc[df["sample"].eq("S14"), "jcpyv_mapq30_first_mates"].iloc[0])),
                    xytext=(11, 180), arrowprops={"arrowstyle": "->", "color": "#555555"})
    axs[2].plot(x, df.bkpyv_breadth_100x_pct, "o", color="#236c91", label="BKPyV ≥100×")
    axs[2].plot(x, df.jcpyv_breadth_100x_pct, "s", color="#a94a34", label="JCPyV ≥100×")
    axs[2].set(ylim=(-4, 105), ylabel="Genome breadth (%)")
    axs[2].legend(ncol=2, frameon=False, loc="center left")
    axs[3].plot(x, df.shared_bk_jc_3p_mirna_raw_reads + 1, "o", color="#7650a0",
                label="Shared BK/JC 3p")
    axs[3].plot(x, df.bkpyv_5p_mirna_raw_reads + 1, "s", color="#25836c",
                label="BKPyV 5p")
    axs[3].set(yscale="log", ylabel="Small-RNA raw\nreads + 1")
    axs[3].legend(ncol=2, frameon=False, loc="upper left")
    axs[-1].set(xticks=x, xticklabels=df["sample"], xlabel="Blood-DNAemia-negative controls (left)  |  cases (right)")
    fig.suptitle("Urinary polyomavirus DNA evidence and separately scaled viral miRNAs", fontsize=12)
    save(fig, figdir, "Figure_2_polyomavirus")

    fig, axs = plt.subplots(3, 1, figsize=(12, 6.7), sharex=True, layout="constrained")
    for ax in axs:
        ax.axvspan(10.5, 21.5, color="#a94a34", alpha=.035)
        ax.axvline(10.5, color="#777777", linewidth=.8)
        ax.grid(axis="y", color="#dddddd", linewidth=.5)
    axs[0].scatter(x, df.smallrna_retained_reads + 1, c=colors, s=30)
    axs[0].set(yscale="log", ylabel="Retained small-RNA\nreads + 1")
    axs[1].scatter(x, df.host_filtered_dna_pairs + 1, c=colors, s=30)
    axs[1].set(yscale="log", ylabel="Host-filtered DNA\nread pairs + 1")
    axs[2].scatter(x, df.kraken_human_pct_classifier_input, c=colors, s=30)
    axs[2].set_yscale("symlog", linthresh=.01)
    axs[2].set_ylabel("Residual human Kraken2\nassignments (% input)")
    axs[-1].set(xticks=x, xticklabels=df["sample"],
                xlabel="Blood-DNAemia-negative controls (left)  |  cases (right)")
    fig.suptitle("Sample yields and post-host-filter human-read screen", fontsize=12)
    save(fig, figdir, "Supplement_sample_qc")

    profiles = a.root / "results/tables"
    sample_groups = ma.loc[order, "group"]
    ranks = [("class", ["Gammaproteobacteria", "Actinomycetes", "Bacilli"]),
             ("genus", ["Serratia", "Corynebacterium", "Streptococcus"]),
             ("species", ["Serratia_marcescens", "Escherichia_coli", "Lactobacillus_iners"])]
    fig, axs = plt.subplots(3, 3, figsize=(12, 8), sharex=True, layout="constrained")
    bacterial_rows = []
    for row, (rank, taxa) in enumerate(ranks):
        profile = pd.read_csv(profiles / f"{rank}_profiles_all22.csv", index_col=0)
        assert set(order) == set(profile.columns)
        for col, taxon in enumerate(taxa):
            vals = profile.loc[taxon, order].astype(float)
            bacterial_rows.extend({"rank": rank, "taxon": taxon, "sample": sample,
                                   "group": sample_groups[sample], "relative_pct": vals[sample]}
                                  for sample in order)
            ax = axs[row, col]
            ax.scatter(x, vals, c=colors, s=26, zorder=3)
            ax.axvline(10.5, color="#777777", linewidth=.8)
            ax.axvspan(10.5, 21.5, color="#a94a34", alpha=.035)
            ax.set_title(taxon.replace("_", " "), fontsize=9)
            ax.set_ylim(bottom=-.8)
            ax.grid(axis="y", color="#dddddd", linewidth=.5)
            if col == 0:
                ax.set_ylabel(f"{rank.title()}\nbacterial relative %")
            if row == 2:
                ax.set_xticks(x)
                ax.set_xticklabels(df["sample"], rotation=90)
    fig.suptitle("Individual bacterial profiles: classes and representative lower taxa", fontsize=12)
    save(fig, figdir, "Figure_3_bacterial_individuals")
    pd.DataFrame(bacterial_rows).to_csv(tabdir / "plotted_bacterial_sample_values.csv", index=False)
    summary = {
        "n_samples": len(df), "n_cases": int(df.group.eq("case").sum()),
        "n_controls": int(df.group.eq("control").sum()),
        "bkpyv_100x_breadth_min_pct": float(df.bkpyv_breadth_100x_pct.min()),
        "bkpyv_high_mapq_first_mates_min": int(df.bkpyv_mapq30_first_mates.min()),
        "jcpyv_100x_samples": df.loc[df.jcpyv_breadth_100x_pct > 0, "sample"].tolist(),
        "jcpyv_max_other_10x_breadth_pct": float(df.loc[df["sample"].ne("S14"), "jcpyv_breadth_10x_pct"].max()),
        "figure_note": "Descriptive DNA/RNA integration after reviewing polyomavirus alignment; not a prespecified association test. Blood DNAemia defines the groups. DNA and small-RNA counts have separate denominators.",
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
