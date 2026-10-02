#!/usr/bin/env python3
"""Assemble Figure 1 differential miRNAs, Figure 2 associations, Figure 3 Sankey."""
import argparse
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from sankey_figure import draw_sankey


HISTORICAL = [
    "hsa-miR-345-5p", "hsa-miR-361-5p", "hsa-miR-193a-3p", "hsa-miR-339-3p",
    "hsa-miR-181c-3p", "hsa-miR-96-5p", "hsa-miR-151a-5p", "hsa-miR-320b",
    "hsa-miR-589-5p", "hsa-let-7a-3p", "bkv-miR-B1-3p", "bkv-miR-B1-5p",
]


def copy_figure(source, target, new_stem=None):
    stem = source.stem if new_stem is None else new_stem
    for extension in ("png","svg"):
        shutil.copy2(source.with_suffix("."+extension),target/(stem+"."+extension))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root","joint-figure","heatmap","previous-figures","out"):
        parser.add_argument("--"+name,type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"Output directory exists: {args.out}")
    args.out.mkdir(parents=True)
    copy_figure(args.joint_figure / "Figure_1_joint_human_viral_miRNA",args.out)
    copy_figure(args.heatmap / "Figure_2_miRNA_genus_heatmap",args.out)
    for n in (1,2):
        source = args.joint_figure if n==1 else args.heatmap
        shutil.copy2(source/f"Figure_{n}_caption.txt",args.out/f"Figure_{n}_caption.txt")
    copy_figure(args.previous_figures / "Figure_2_polyomavirus",args.out,
                "Supplement_polyomavirus_DNA_RNA")
    copy_figure(args.previous_figures / "Figure_3_bacterial_individuals",args.out,
                "Supplement_bacterial_individuals")
    copy_figure(args.previous_figures / "Supplement_historical_human_miRNAs",args.out)

    metadata = pd.read_csv(args.root/"verified_data/metadata/bkv_sample_metadata.csv").set_index("sample")
    sf = pd.read_csv(args.root/"results/count_models/main_ge2000/size_factors.csv").set_index("sample").size_factor
    samples = sf.index.tolist()
    controls = [s for s in samples if metadata.loc[s,"group"]=="control"]
    cases = [s for s in samples if metadata.loc[s,"group"]=="treatment"]
    assert len(controls)==7 and len(cases)==9
    classes = pd.read_csv(args.root/"results/tables/class_profiles_all22.csv",index_col=0)
    human = pd.read_csv(args.root/"inputs/recount/human_unique_counts.csv",index_col=0)
    viral = pd.read_csv(args.root/"inputs/recount/viral_historical_counts.csv",index_col=0)
    normalized = np.log2(pd.concat([human,viral])[samples].div(sf)+1).loc[HISTORICAL]
    draw_sankey(classes[samples],normalized,controls,cases,args.out,args.out,
                figure_number=3)
    (args.out/"Figure_3_caption.txt").write_text(
        "Figure 3. Descriptive Sankey for the paired 9-case/7-control subset. "
        "Left ribbons encode bacterial class percentages; classes remaining below 1% "
        "in every displayed sample are grouped as Other. Dashed guides make small Other "
        "connections visible but do not encode abundance. Right ribbons encode absolute "
        "deviations of log2(count/human DESeq2 size factor+1) from the control mean, "
        "displayed when the deviation is at least 0.5. Human counts use revised "
        "unique-mature assignments; viral counts retain BK-only attribution. Controls "
        "contribute to their reference mean. Human above/below deviations are red/green, "
        "viral above/below are dark/light blue. Each side has an independent linear "
        "width scale. Links indicate sample membership, not correlations or interactions. "
        "The asterisk marks the BKPyV/JCPyV-shared 3p sequence.\n"
    )
    print(args.out)


if __name__ == "__main__":
    main()
