#!/usr/bin/env python3
"""Summarize saved Kraken2 reports without re-reading large FASTQ/assignment files.

The input is the small_outputs.tar.gz archive written by run_dna_kraken2.sbatch
for each accession. Counts refer to paired-read classifications, not organism
abundance per mL, clinical infection, or validated species-level identifications.
"""

import argparse
import csv
import hashlib
import re
import tarfile
from pathlib import Path


TAXA = {
    "bacteria": "2",
    "archaea": "2157",
    "viruses": "10239",
    "eukaryota": "2759",
    "fungi": "4751",
    "human": "9606",
    "bkpyv": "1891762",
    "jcpyv": "10632",
}


def read_table(path, delimiter):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_archive(path):
    with tarfile.open(path, "r:gz") as archive:
        files = {}
        for name in ("report.tsv", "run_metadata.txt", "kraken_stderr.txt", "SUCCESS"):
            member = archive.getmember(name)
            if not member.isfile():
                raise ValueError(f"{path}: {name} is not a regular file")
            files[name] = archive.extractfile(member).read().decode("utf-8")
    if "Complete:" not in files["SUCCESS"]:
        raise ValueError(f"{path}: SUCCESS marker is malformed")
    return files


def metadata_fields(text):
    result = {}
    for line in text.splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            result[key] = value
    return result


def parse_report(text):
    rows, ancestors = [], []
    for line in text.splitlines():
        parts = line.split("\t", 7)
        if len(parts) != 8:
            raise ValueError(f"Expected eight Kraken2 report columns: {line[:100]!r}")
        percent, clade, direct, minimizers, unique, rank, taxid, indented_name = parts
        depth = len(indented_name) - len(indented_name.lstrip(" "))
        while ancestors and ancestors[-1][0] >= depth:
            ancestors.pop()
        name = indented_name.strip()
        lineage = [item[1] for item in ancestors] + [name]
        ancestors.append((depth, name))
        rows.append({
            "percent_reported": percent.strip(),
            "clade_pairs": int(clade),
            "direct_pairs": int(direct),
            "minimizers": int(minimizers),
            "distinct_minimizers": int(unique),
            "rank": rank,
            "taxid": taxid,
            "name": name,
            "lineage": "; ".join(lineage),
        })
    return rows


def parse_stderr(text):
    match = re.search(r"^(\d+) sequences \(.+\) processed in ([0-9.]+)s", text, re.M)
    classified = re.search(r"^\s*(\d+) sequences classified", text, re.M)
    unclassified = re.search(r"^\s*(\d+) sequences unclassified", text, re.M)
    if not (match and classified and unclassified):
        raise ValueError("Cannot parse Kraken2 total/classified/unclassified read-pair counts")
    total, seconds = int(match.group(1)), float(match.group(2))
    c, u = int(classified.group(1)), int(unclassified.group(1))
    if c + u != total:
        raise ValueError("Kraken2 classified + unclassified does not equal total")
    return total, c, u, seconds


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-dir", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=22)
    args = parser.parse_args()
    if args.out_dir.exists():
        parser.error(f"Output directory already exists: {args.out_dir}")

    mapping_rows = read_table(args.mapping, "\t")
    accession_to_sample = {r["run_accession"]: f'S{r["sample_title"]}' for r in mapping_rows}
    if len(accession_to_sample) != len(mapping_rows):
        raise ValueError("Duplicate accessions in mapping")
    group_rows = read_table(args.metadata, ",")
    sample_to_group = {r["sample"]: r["group"] for r in group_rows}
    if len(sample_to_group) != len(group_rows):
        raise ValueError("Duplicate samples in clinical metadata")

    archives = sorted(args.reports_dir.glob("SRR*_small_outputs.tar.gz"))
    if len(archives) != args.expected:
        raise ValueError(f"Expected {args.expected} reports, found {len(archives)} in {args.reports_dir}")
    summary_rows, taxon_rows = [], []
    seen, db_signatures = set(), set()
    for path in archives:
        accession = path.name.split("_small_outputs.tar.gz")[0]
        if accession in seen or accession not in accession_to_sample:
            raise ValueError(f"Duplicate or unmapped accession: {accession}")
        seen.add(accession)
        files = read_archive(path)
        meta = metadata_fields(files["run_metadata.txt"])
        if meta.get("sample") != accession:
            raise ValueError(f"Sample mismatch in {path}: {meta.get('sample')}")
        if meta.get("memory_mapping") != "0":
            raise ValueError(f"Unexpected memory mapping setting in {path}")
        if "Kraken version 2.1.2" not in files["run_metadata.txt"]:
            raise ValueError(f"Unexpected classifier version in {path}")
        db_signatures.add((meta.get("database"), files["run_metadata.txt"].split("database_files\n", 1)[1].split("input_files\n", 1)[0]))
        total, classified, unclassified, seconds = parse_stderr(files["kraken_stderr.txt"])
        report_rows = parse_report(files["report.tsv"])
        by_taxid = {r["taxid"]: r for r in report_rows}
        if len(by_taxid) != len(report_rows):
            raise ValueError(f"Duplicate taxids in report for {accession}")
        if by_taxid["0"]["clade_pairs"] != unclassified or by_taxid["1"]["clade_pairs"] != classified:
            raise ValueError(f"Report totals disagree with stderr for {accession}")
        sample = accession_to_sample[accession]
        if sample not in sample_to_group:
            raise ValueError(f"Missing group metadata for {sample}")
        row = {
            "accession": accession, "sample": sample, "group": sample_to_group[sample],
            "total_pairs": total, "classified_pairs": classified,
            "unclassified_pairs": unclassified, "classifier_seconds": seconds,
            "report_sha256": hashlib.sha256(files["report.tsv"].encode()).hexdigest(),
            "database": meta.get("database", ""), "kraken_version": "2.1.2",
        }
        for label, taxid in TAXA.items():
            taxon = by_taxid.get(taxid)
            row[f"{label}_pairs"] = taxon["clade_pairs"] if taxon else 0
            row[f"{label}_distinct_minimizers"] = taxon["distinct_minimizers"] if taxon else 0
        summary_rows.append(row)
        for taxon in report_rows:
            if taxon["rank"] == "U" or taxon["clade_pairs"] == 0:
                continue
            taxon_rows.append({"accession": accession, "sample": sample, "group": sample_to_group[sample], **taxon})

    if len(db_signatures) != 1:
        raise ValueError("Reports used differing database paths or database file signatures")
    args.out_dir.mkdir(parents=True)
    write_tsv(args.out_dir / "kraken2_sample_summary.tsv", summary_rows, list(summary_rows[0]))
    write_tsv(args.out_dir / "kraken2_all_taxa.tsv", taxon_rows, list(taxon_rows[0]))
    print(f"Validated {len(archives)} reports; wrote {len(taxon_rows)} nonzero taxon rows to {args.out_dir}")


if __name__ == "__main__":
    main()
