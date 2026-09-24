#!/usr/bin/env python3
"""Export descriptive per-sample fungal and other eukaryotic Kraken2 labels."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxa", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        p.error(f"Output already exists: {a.out}")
    a.out.mkdir(parents=True)
    x = pd.read_csv(a.taxa, sep="\t")
    s = pd.read_csv(a.summary, sep="\t")
    assert len(s) == 22 and s.accession.is_unique and s["sample"].is_unique
    e = x[x["rank"].eq("S") & x.lineage.str.contains("Eukaryota", regex=False)].copy()
    e["screen_category"] = "other_eukaryote"
    e.loc[e.lineage.str.contains("Fungi", regex=False), "screen_category"] = "fungus"
    e.loc[e.lineage.str.contains("Metazoa", regex=False), "screen_category"] = "metazoan"
    e.loc[e.lineage.str.contains("Viridiplantae", regex=False), "screen_category"] = "green_plant"
    e = e[e.screen_category.isin(["fungus", "other_eukaryote"])]
    e = e.merge(s[["accession", "sample", "group", "total_pairs"]],
                on=["accession", "sample", "group"], validate="many_to_one")
    e["pct_classifier_input"] = e.clade_pairs / e.total_pairs * 100
    cols = ["accession", "sample", "group", "screen_category", "taxid", "name",
            "clade_pairs", "direct_pairs", "distinct_minimizers", "total_pairs",
            "pct_classifier_input", "lineage"]
    e[cols].sort_values(["sample", "screen_category", "clade_pairs"],
                        ascending=[True, True, False]).to_csv(a.out / "eukaryotic_species_labels.tsv",
                                                              sep="\t", index=False)
    agg = e.groupby(["sample", "screen_category"]).clade_pairs.sum().unstack(fill_value=0)
    agg = s.set_index("sample")[["accession", "group", "total_pairs", "fungi_pairs",
                                 "fungi_distinct_minimizers", "human_pairs"]].join(agg)
    for category in ["fungus", "other_eukaryote"]:
        if category not in agg:
            agg[category] = 0
        agg[category + "_species_clade_pairs_pct_input"] = agg[category] / agg.total_pairs * 100
    agg["kraken_fungi_clade_pct_input"] = agg.fungi_pairs / agg.total_pairs * 100
    agg.sort_index().to_csv(a.out / "eukaryote_screen_by_sample.tsv", sep="\t")
    info = {
        "taxa_sha256": hashlib.sha256(a.taxa.read_bytes()).hexdigest(),
        "summary_sha256": hashlib.sha256(a.summary.read_bytes()).hexdigest(),
        "interpretation": "Kraken2 PlusPF taxon labels only. Other eukaryote denotes Eukaryota species lineage without Fungi, Metazoa or Viridiplantae. Summed species clade pairs are a descriptive screen, not independent infection evidence or exhaustive parasite detection."
    }
    (a.out / "provenance.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"Exported {len(e)} species-label rows across {e['sample'].nunique()} samples")


if __name__ == "__main__":
    main()
