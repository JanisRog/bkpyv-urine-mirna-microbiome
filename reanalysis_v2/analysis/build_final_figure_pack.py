#!/usr/bin/env python3
"""Assemble numbered main figures while retaining the descriptive Sankey."""

import argparse
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
from sankey_figure import draw_sankey


p = argparse.ArgumentParser(description=__doc__)
for key in ["root", "table-run", "dna-integration", "out"]:
    p.add_argument("--" + key, type=Path, required=True)
a = p.parse_args()
if a.out.exists():
    p.error(f"Output exists: {a.out}")
a.out.mkdir(parents=True)
for src, stem in [
    (a.table_run, "Figure_1_human_miRNA"),
    (a.dna_integration, "Figure_2_polyomavirus"),
    (a.dna_integration, "Figure_3_bacterial_individuals"),
]:
    for suffix in ["png", "svg"]:
        shutil.copy2(src / "figures" / f"{stem}.{suffix}", a.out / f"{stem}.{suffix}")

meta = pd.read_csv(a.root / "verified_data/metadata/bkv_sample_metadata.csv").set_index("sample")
sf = pd.read_csv(a.table_run / "results/count_models/main_ge2000/size_factors.csv").set_index("sample").size_factor
ss = sf.index.tolist()
controls = [s for s in ss if meta.loc[s, "group"] == "control"]
cases = [s for s in ss if meta.loc[s, "group"] == "treatment"]
cls = pd.read_csv(a.table_run / "results/tables/class_profiles_all22.csv", index_col=0)
human = pd.read_csv(a.root / "inputs/recount/human_unique_counts.csv", index_col=0)
viral = pd.read_csv(a.root / "inputs/recount/viral_historical_counts.csv", index_col=0)
names = ["hsa-miR-345-5p", "hsa-miR-361-5p", "hsa-miR-193a-3p", "hsa-miR-339-3p",
         "hsa-miR-181c-3p", "hsa-miR-96-5p", "hsa-miR-151a-5p", "hsa-miR-320b",
         "hsa-miR-589-5p", "hsa-let-7a-3p", "bkv-miR-B1-3p", "bkv-miR-B1-5p"]
normalized = np.log2(pd.concat([human, viral])[ss].div(sf) + 1).loc[names]
draw_sankey(cls[ss], normalized, controls, cases, a.out, a.out, figure_number=4)
print("Final Figure 1–4 pack:", a.out)
