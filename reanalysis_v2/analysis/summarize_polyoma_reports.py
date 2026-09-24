#!/usr/bin/env python3
"""Validate per-sample BK/JC alignment archives and assemble a study-wide table."""

import argparse
import csv
import hashlib
import io
import statistics
import tarfile
from pathlib import Path


REFS = {"NC_001538.1": (5153, "bkpyv"), "NC_001699.1": (5130, "jcpyv")}
REF_SHA256 = "26dc1b8c7117a6e52b72cc17dba547887f7e660b4fb2aeae7e38076fdfafdb23"


def table(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def archive_text(archive, name):
    member = archive.getmember(name)
    if not member.isfile():
        raise ValueError(f"{name} is not a regular archive member")
    return archive.extractfile(member).read().decode("utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reports-dir", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--kraken-summary", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--expected", type=int, default=22)
    a = p.parse_args()
    if a.out.exists():
        p.error(f"Output already exists: {a.out}")
    mapping = {x["run_accession"]: "S" + x["sample_title"] for x in table(a.mapping)}
    groups = {x["sample"]: x["group"] for x in table(a.metadata, ",")}
    kraken = {x["accession"]: x for x in table(a.kraken_summary)}
    paths = sorted(a.reports_dir.glob("SRR*_polyoma_summary.tar.gz"))
    if len(paths) != a.expected:
        raise ValueError(f"Expected {a.expected} archives; found {len(paths)}")
    rows, jobs, seen = [], set(), set()
    for path in paths:
        accession = path.name.split("_polyoma_summary.tar.gz")[0]
        if accession in seen or accession not in mapping or accession not in kraken:
            raise ValueError(f"Duplicate or unmapped accession: {accession}")
        seen.add(accession)
        with tarfile.open(path, "r:gz") as archive:
            success = archive_text(archive, "SUCCESS")
            metadata = archive_text(archive, "run_metadata.txt")
            summary = archive_text(archive, "summary.tsv")
            depth = archive_text(archive, "depth_mapq30_baseq20.tsv")
        if success.strip() != f"Complete: {accession}" or f"sample\t{accession}" not in metadata:
            raise ValueError(f"Unsuccessful or mislabeled job: {accession}")
        fields = dict(line.split("\t", 1) for line in metadata.splitlines() if "\t" in line)
        jobs.add(fields["array_job_id"])
        if fields.get("reference_sha256") != REF_SHA256:
            raise ValueError(f"Reference mismatch for {accession}")
        sample = mapping[accession]
        if sample not in groups:
            raise ValueError(f"Missing group for {sample}")
        summary_rows = list(csv.DictReader(io.StringIO(summary), delimiter="\t"))
        if {x["reference"] for x in summary_rows} != set(REFS):
            raise ValueError(f"Unexpected reference list in {accession}")
        depths = {ref: [] for ref in REFS}
        for line in depth.splitlines():
            ref, pos, value = line.split("\t")
            if ref not in depths or int(pos) != len(depths[ref]) + 1:
                raise ValueError(f"Unsorted/invalid depth line for {accession}: {line[:60]}")
            depths[ref].append(int(value))
        for x in summary_rows:
            ref = x["reference"]
            length, label = REFS[ref]
            values = depths[ref]
            if len(values) != length or int(x["reference_bp"]) != length or x["reference_sha256"] != REF_SHA256:
                raise ValueError(f"Depth/reference mismatch: {accession} {ref}")
            for threshold in (1, 10, 100):
                field = f"breadth_{threshold}x_pct_mapq30_baseq20"
                observed = round(100 * sum(v >= threshold for v in values) / length, 3)
                if abs(observed - float(x[field])) > 0.001:
                    raise ValueError(f"Coverage mismatch: {accession} {ref} {field}")
            if abs(statistics.median(values) - float(x["median_depth_mapq30_baseq20"])) > 0.001:
                raise ValueError(f"Median depth mismatch: {accession} {ref}")
            if int(x["mapq30_first_mates"]) > int(x["primary_mapped_alignments"]):
                raise ValueError(f"High-MAPQ count exceeds mapped alignments: {accession} {ref}")
            rows.append({
                "accession": accession, "sample": sample, "group": groups[sample],
                "total_host_filtered_pairs": kraken[accession]["total_pairs"],
                "kraken_species_pairs": kraken[accession][f"{label}_pairs"],
                "kraken_distinct_minimizers": kraken[accession][f"{label}_distinct_minimizers"],
                **{key: value for key, value in x.items() if key != "sample"},
                "alignment_archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "job_id": fields["array_job_id"],
            })
    if len(jobs) != 1:
        raise ValueError(f"Multiple array job IDs: {jobs}")
    with a.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Validated {len(paths)} archives from job {next(iter(jobs))}; wrote {len(rows)} reference rows to {a.out}")


if __name__ == "__main__":
    main()
