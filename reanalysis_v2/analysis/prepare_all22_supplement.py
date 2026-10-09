"""Prepare exact all-sample matrices for a supplementary workbook."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "paired_primary_2026-10-09"
DEST.mkdir(exist_ok=True)

meta = pd.read_csv(ROOT / "verified_data/metadata/bkv_sample_metadata.csv")
samples = sorted(meta["sample"], key=lambda s: int(s[1:]))
assert len(samples) == 22 and len(set(samples)) == 22
primary = set(pd.read_csv(
    ROOT / "joint_mirna_models_2026-10-02_v2/main_ge2000/fixed_host_normalization.csv"
)["sample"])
assert len(primary) == 16 and primary.issubset(samples)
rna_qc = pd.read_csv(ROOT / "inputs/smallRNA_sample_qc.csv").set_index("sample")
human = pd.read_csv(ROOT / "inputs/recount/human_unique_counts.csv", index_col=0)[samples]
viral = pd.read_csv(ROOT / "inputs/recount/viral_historical_counts.csv", index_col=0)[samples]
mapping = pd.read_csv(ROOT / "results/dna_qc/ena_sample_mapping.tsv", sep="\t")
mapping["sample"] = "S" + mapping["sample_title"].astype(str)
dna_qc = pd.read_csv(ROOT / "results/dna_qc/dna_read_qc_by_accession.csv")
dna = mapping.merge(dna_qc, left_on="run_accession", right_on="accession", validate="one_to_one").set_index("sample")
assert set(dna.index) == set(samples)
assert dna["processed_matches_both_unmapped"].all()
assert ((dna["species_sum"] - 100).abs() < 0.001).all()
assert human.sum(axis=0).ge(2000).to_dict() == {s: s in primary for s in samples}

sample_rows = [[
    "Sample ID", "Blood BKPyV DNAemia group", "Primary paired cohort",
    "Human mature miRNA assigned counts", "R1 input reads", "R1 retained reads",
    "MetaPhlAn input reads", "Primary human DNA aligned (%)",
    "Species in MetaPhlAn profile",
]]
for sample in samples:
    group = meta.set_index("sample").loc[sample, "group"]
    sample_rows.append([
        sample, "Case" if group == "treatment" else "Control",
        "Yes" if sample in primary else "No",
        int(human[sample].sum()), int(rna_qc.loc[sample, "input_reads"]),
        int(rna_qc.loc[sample, "retained_reads"]),
        int(dna.loc[sample, "metaphlan_reported_reads"]),
        float(dna.loc[sample, "host_mapped_pct"]),
        int(dna.loc[sample, "profile_species"]),
    ])


def matrix_rows(frame: pd.DataFrame, first_heading: str) -> list[list[object]]:
    assert frame.columns.tolist() == samples
    rows = [[first_heading, *samples]]
    for feature, values in frame.iterrows():
        rows.append([str(feature), *[float(x) if isinstance(x, float) else int(x) for x in values]])
    return rows


matrices = {
    "Human miRNA counts": matrix_rows(human, "Mature miRNA"),
    "Viral miRNA counts": matrix_rows(viral, "Viral mature miRNA"),
}
for rank in ("class", "genus", "species"):
    frame = pd.read_csv(ROOT / f"results/tables/{rank}_profiles_all22.csv", index_col=0)[samples]
    assert ((frame.sum(axis=0) - 100).abs() < 0.001).all()
    matrices[f"Bacterial {rank} %"] = matrix_rows(frame, rank.title())

payload = {
    "sampleRows": sample_rows,
    "matrices": matrices,
    "notes": [
        ["Supplementary Table S6", "Urinary host and viral mature miRNA counts and MetaPhlAn bacterial relative percentages for all 22 recipients."],
        ["Primary cohort", "Nine blood-DNAemia cases and seven controls with at least 2,000 uniquely assigned mature human miRNA counts. This retrospective yield threshold applies to the paired primary analysis."],
        ["Six outside primary cohort", "S1, S11, S12, S13, S15 and S20 are retained in the all-22 sensitivity data. Their low miRNA yield is not a failure of the DNA/MetaPhlAn profile checks."],
        ["Counts", "Human and viral miRNA sheets contain per-sample integer mature-sequence counts before between-sample normalization. Zero means no assigned reads."],
        ["Shared viral sequence", "bkv-miR-B1-3p is identical in sequence to jcv-miR-J1-3p and cannot distinguish the two viruses."],
        ["Bacterial percentages", "Class, genus and species values are MetaPhlAn relative percentages among bacterial profile assignments. Each sample sums to approximately 100% within each rank. These are not read counts or organisms per mL."],
        ["Sample codes", "S-prefixed identifiers are deidentified analysis sample codes; cases had blood BKPyV DNAemia, controls did not."],
        ["Analysis code", "https://github.com/JanisRog/bkpyv-urine-mirna-microbiome"],
    ],
}
(DEST / "all22_supplement_data.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
print(json.dumps({"samples": len(samples), "human_features": len(human), "viral_features": len(viral), "bacterial_features": {k: len(v)-1 for k,v in matrices.items() if k.startswith("Bacterial")}}))
