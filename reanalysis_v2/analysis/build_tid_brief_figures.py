#!/usr/bin/env python3
"""Assemble the two main figures and three supplementary figures for TID."""
import argparse
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from sankey_figure import draw_sankey


HISTORICAL_MIRNAS = [
    "hsa-miR-345-5p", "hsa-miR-361-5p", "hsa-miR-193a-3p",
    "hsa-miR-339-3p", "hsa-miR-181c-3p", "hsa-miR-96-5p",
    "hsa-miR-151a-5p", "hsa-miR-320b", "hsa-miR-589-5p",
    "hsa-let-7a-3p", "bkv-miR-B1-3p", "bkv-miR-B1-5p",
]

VIRAL_COLUMNS = [
    "sample", "group", "host_filtered_dna_pairs", "bkpyv_mapq30_first_mates",
    "jcpyv_mapq30_first_mates", "bkpyv_mapq30_first_mates_pct_host_filtered_pairs",
    "jcpyv_mapq30_first_mates_pct_host_filtered_pairs", "bkpyv_breadth_100x_pct",
    "jcpyv_breadth_100x_pct", "shared_bk_jc_3p_mirna_raw_reads",
    "bkpyv_5p_mirna_raw_reads",
]


def copy_figure(source, directory, stem):
    for extension in ("png", "svg"):
        shutil.copy2(source.with_suffix("." + extension), directory / (stem + "." + extension))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "joint-figure", "heatmap", "integration", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"Output directory already exists: {args.out}")
    main_dir = args.out / "main"
    supplement = args.out / "supplementary"
    main_dir.mkdir(parents=True)
    supplement.mkdir()

    copy_figure(args.joint_figure / "Figure_1_joint_human_viral_miRNA", main_dir,
                "Figure_1_joint_human_viral_miRNA")
    copy_figure(args.heatmap / "Figure_2_miRNA_genus_heatmap", main_dir,
                "Figure_2_miRNA_genus_heatmap")
    copy_figure(args.integration / "figures/Figure_2_polyomavirus", supplement,
                "Figure_S1_polyomavirus_DNA_RNA")
    copy_figure(args.integration / "figures/Figure_3_bacterial_individuals", supplement,
                "Figure_S2_bacterial_individuals")

    meta = pd.read_csv(args.root / "verified_data/metadata/bkv_sample_metadata.csv").set_index("sample")
    factors = pd.read_csv(args.root / "results/count_models/main_ge2000/size_factors.csv").set_index("sample").size_factor
    samples = factors.index.tolist()
    controls = [sample for sample in samples if meta.loc[sample, "group"] == "control"]
    cases = [sample for sample in samples if meta.loc[sample, "group"] == "treatment"]
    assert len(controls) == 7 and len(cases) == 9
    classes = pd.read_csv(args.root / "results/tables/class_profiles_all22.csv", index_col=0)
    human = pd.read_csv(args.root / "inputs/recount/human_unique_counts.csv", index_col=0)
    viral = pd.read_csv(args.root / "inputs/recount/viral_historical_counts.csv", index_col=0)
    normalized = np.log2(pd.concat([human, viral])[samples].div(factors) + 1).loc[HISTORICAL_MIRNAS]
    draw_sankey(classes[samples], normalized, controls, cases, supplement, supplement,
                figure_number="S3")

    viral_data = pd.read_csv(args.integration / "tables/sample_level_viral_and_qc.csv")
    assert len(viral_data) == 22 and set(VIRAL_COLUMNS) <= set(viral_data.columns)
    viral_data[VIRAL_COLUMNS].to_csv(
        supplement / "Supplementary_Table_S4_polyomavirus_reads.csv", index=False)
    shutil.copy2(args.heatmap / "Figure_2_miRNA_genus_association_tests.csv",
                 supplement / "Supplementary_Table_S5_miRNA_genus_associations.csv")
    print(args.out)


if __name__ == "__main__":
    main()
