#!/usr/bin/env python3
"""Validate and summarize BK/JC competitive-alignment depth and read counts."""

import argparse
import csv
import hashlib
import statistics
from pathlib import Path


EXPECTED = {"NC_001538.1": 5153, "NC_001699.1": 5130}


def fasta_lengths(path):
    lengths = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in lengths:
                raise ValueError(f"Duplicate FASTA record: {name}")
            lengths[name] = 0
        elif name:
            lengths[name] += len(line.strip())
    if lengths != EXPECTED:
        raise ValueError(f"Unexpected reference lengths: {lengths}")
    return lengths


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--depth", type=Path, required=True)
    p.add_argument("--idxstats", type=Path, required=True)
    p.add_argument("--counts", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()

    lengths = fasta_lengths(a.reference)
    depths = {ref: [None] * length for ref, length in lengths.items()}
    with a.depth.open() as handle:
        for line in handle:
            ref, pos, value = line.rstrip("\n").split("\t")
            if ref not in depths:
                raise ValueError(f"Unexpected depth reference: {ref}")
            pos = int(pos)
            if pos < 1 or pos > lengths[ref] or depths[ref][pos - 1] is not None:
                raise ValueError(f"Invalid/duplicate depth position: {ref}:{pos}")
            depths[ref][pos - 1] = int(value)
    if any(None in values for values in depths.values()):
        raise ValueError("Depth table omits reference positions; use samtools depth -aa")

    idxstats = {}
    with a.idxstats.open() as handle:
        for line in handle:
            ref, length, mapped, unmapped = line.rstrip("\n").split("\t")
            if ref == "*":
                continue
            idxstats[ref] = (int(length), int(mapped), int(unmapped))
    if set(idxstats) != set(lengths) or any(idxstats[r][0] != lengths[r] for r in lengths):
        raise ValueError(f"Unexpected idxstats reference names/lengths: {idxstats}")

    counts = {}
    with a.counts.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            counts[row["reference"]] = int(row["mapq30_first_mates"])
    if set(counts) != set(lengths):
        raise ValueError(f"Missing high-MAPQ counts for: {set(lengths) - set(counts)}")

    rows = []
    for ref, length in lengths.items():
        values = depths[ref]
        rows.append({
            "sample": a.sample,
            "reference": ref,
            "reference_bp": length,
            "reference_sha256": hashlib.sha256(a.reference.read_bytes()).hexdigest(),
            "primary_mapped_alignments": idxstats[ref][1],
            "mapq30_first_mates": counts[ref],
            "breadth_1x_pct_mapq30_baseq20": round(100 * sum(v >= 1 for v in values) / length, 3),
            "breadth_10x_pct_mapq30_baseq20": round(100 * sum(v >= 10 for v in values) / length, 3),
            "breadth_100x_pct_mapq30_baseq20": round(100 * sum(v >= 100 for v in values) / length, 3),
            "median_depth_mapq30_baseq20": statistics.median(values),
            "mean_depth_mapq30_baseq20": round(statistics.mean(values), 3),
            "max_depth_mapq30_baseq20": max(values),
        })
    with a.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
