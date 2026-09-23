#!/usr/bin/env python3
"""Recount existing miRDeep2 mappings without repeating read alignment.

Standard library only. Preserves original read multiplicities. Collapses repeat
assignments to the SAME mature name; ambiguous DIFFERENT mature names are
excluded from unique-mature counts and retained in explicit ambiguity groups.
Historical precursor counts must reconstruct exactly before any matrix is saved.
"""
import argparse
from collections import Counter, defaultdict
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import re
import sys


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(path, header, rows):
    with path.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def integer(value):
    n = Decimal(value)
    if not n.is_finite() or n < 0 or n != int(n):
        raise ValueError('Expected a finite nonnegative integer count: ' + value)
    return int(n)


def only(paths):
    paths = list(paths)
    if len(paths) != 1:
        raise ValueError('Expected exactly one file: ' + str(paths))
    return paths[0]


def read_expression(path):
    expected = Counter()
    features = set()
    with path.open(newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            mir, pre = row['#miRNA'], row['precursor']
            key = (mir, pre)
            # Multiple mature placements can generate repeated pair rows.
            # Preserve their sum and require exact reconstruction from mappings.
            expected[key] += integer(row['read_count'])
            features.add(mir)
    if not expected:
        raise ValueError('Empty expression table: ' + str(path))
    return expected, features


def arf(path):
    with path.open() as f:
        for line_number, line in enumerate(f, 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 13:
                raise ValueError(f'Unexpected ARF row: {path}:{line_number}')
            if fields[10] != '+':
                raise ValueError(f'Unexpected antisense mapping: {path}:{line_number}')
            yield fields


def recount(expression, mature_arf, reads_arf):
    expected, features = read_expression(expression)
    windows = defaultdict(list)
    for f in arf(mature_arf):
        mir, pre = f[0], f[5]
        if (mir, pre) in expected:
            # ARF is 1-based. Historical miRDeep2 mature allowance: -2/+5 nt.
            windows[pre].append((mir, max(1, int(f[7])-2), int(f[8])+5))
    reconstructed = Counter()
    assignments = defaultdict(set)
    weights = {}
    for f in arf(reads_arf):
        rid, pre = f[0], f[5]
        m = re.fullmatch(r'[A-Za-z]{3}_\d+_x([1-9]\d*)', rid)
        if not m:
            raise ValueError('Unexpected collapsed-read identifier: ' + rid)
        weight = int(m[1])
        start, end = int(f[7]), int(f[8])
        for mir, left, right in windows[pre]:
            if start >= left and end <= right:
                reconstructed[(mir, pre)] += weight
                assignments[rid].add(mir)
                weights[rid] = weight
    differences = [(key, value, reconstructed[key]) for key, value in expected.items()
                   if value != reconstructed[key]]
    if differences:
        raise ValueError('Historical count reconstruction failed; do not use output. '
                         + str(differences[:10]))
    historical, union, unique, fractional, groups = (Counter() for _ in range(5))
    for (mir, pre), value in expected.items():
        historical[mir] += value
    ambiguous_reads = 0
    for rid, names in assignments.items():
        weight = weights[rid]
        names = tuple(sorted(names))
        # A group owns each observed read once, with its original multiplicity.
        groups[names] += weight
        for name in names:
            union[name] += weight
            fractional[name] += weight / len(names)
        if len(names) == 1:
            unique[names[0]] += weight
        else:
            ambiguous_reads += weight
    if sum(unique.values()) + ambiguous_reads != sum(groups.values()):
        raise AssertionError('Read accounting failed')
    metrics = {'verified_precursor_rows': len(expected),
               'historical_count_total': sum(historical.values()),
               'same_mature_union_total': sum(union.values()),
               'unique_mature_total': sum(unique.values()),
               'ambiguous_read_total': ambiguous_reads,
               'distinct_assigned_read_total': sum(groups.values()),
               'same_mature_repeated_assignment_excess':sum(historical.values())-sum(union.values())}
    return features, historical, union, unique, fractional, groups, metrics


def matrix_check(path, values, features, samples):
    with path.open(newline='') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or reader.fieldnames[0] != 'miRNA':
            raise ValueError('Unexpected matrix header: ' + str(path))
        if set(reader.fieldnames[1:]) != set(samples):
            raise ValueError('Matrix sample mismatch: ' + str(path))
        seen = set()
        for row in reader:
            mir = row['miRNA']
            if mir in seen:
                raise ValueError('Duplicate feature in matrix: ' + mir)
            seen.add(mir)
            for sample in samples:
                if integer(row[sample]) != values[sample].get(mir, 0):
                    raise ValueError(f'Matrix mismatch {path}: {mir}, {sample}')
        if seen != features:
            raise ValueError('Matrix feature mismatch: ' + str(path))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--expected-samples', type=int, default=22)
    a = p.parse_args()
    src = a.source.resolve()
    # Refuse reuse; a failed run remains identifiable and cannot masquerade as success.
    a.out.mkdir(parents=True, exist_ok=False)
    checksums = {}
    def record(path):
        checksums[str(path)] = sha256(path)
    folders = [x for x in (src/'mirdeep2').iterdir() if x.is_dir() and re.fullmatch(r'S\d+', x.name)]
    folders.sort(key=lambda x: int(x.name[1:]))
    samples = [x.name for x in folders]
    if len(samples) != a.expected_samples:
        raise ValueError(f'Expected {a.expected_samples} samples; found {samples}')
    all_features = set()
    matrices = {key:{} for key in ('historical','union','unique','fractional')}
    qc, audit, ambiguity, viral = [], [], [], {}
    for folder in folders:
        sample = folder.name
        expression = only(folder.glob('miRNAs_expressed_all_samples*.csv'))
        exp = folder/'expression_analyses'/f'expression_analyses_{sample}'
        paths = [expression, exp/'mature.fa_mapped.arf', exp/f'{sample}.collapsed.fa_mapped.arf']
        for path in paths: record(path)
        features, historical, union, unique, fractional, groups, metrics = recount(*paths)
        all_features.update(features)
        for key, val in zip(matrices, (historical, union, unique, fractional)):
            matrices[key][sample] = val
        qc.append({'sample':sample, **metrics})
        for mir in sorted(features):
            audit.append([sample, mir, historical[mir], union[mir], unique[mir], fractional[mir]])
        for names, count in sorted(groups.items()):
            if len(names)>1: ambiguity.append([sample, '|'.join(names), len(names), count])
        # Verify/copy historical BK-only counts; no claim of BK/JC-specific reassignment.
        bkv_table = only((src/'mirdeep2_bkv'/sample).glob('miRNAs_expressed_all_samples*.csv'))
        record(bkv_table)
        expected, names = read_expression(bkv_table)
        if names != {'bkv-miR-B1-3p','bkv-miR-B1-5p'}:
            raise ValueError('Unexpected BK-only features: '+str(names))
        viral[sample] = Counter()
        for (mir, pre), count in expected.items(): viral[sample][mir] += count
        print(sample + ': historical precursor counts verified', flush=True)
    for path, values, features in [(src/'counts_matrix.csv', matrices['historical'], all_features),
                                   (src/'bkv_counts_matrix.csv', viral, {'bkv-miR-B1-3p','bkv-miR-B1-5p'})]:
        record(path)
        matrix_check(path, values, features, samples)
    for key in matrices:
        write_csv(a.out/f'human_{key}_counts.csv', ['miRNA']+samples,
                  ([mir]+[matrices[key][s][mir] for s in samples] for mir in sorted(all_features)))
    write_csv(a.out/'viral_historical_counts.csv', ['miRNA']+samples,
              ([mir]+[viral[s][mir] for s in samples] for mir in sorted(viral[samples[0]])))
    write_csv(a.out/'sample_qc.csv', list(qc[0]), ([row[k] for k in qc[0]] for row in qc))
    write_csv(a.out/'feature_audit.csv', ['sample','miRNA','historical','union','unique','fractional'], audit)
    write_csv(a.out/'ambiguous_mature_groups.csv', ['sample','mature_names','number_of_names','read_count'], ambiguity)
    record(Path(__file__).resolve())
    for directory in ('ref','ref_bkv'):
        for path in sorted((src/directory).glob('*.fa')): record(path)
    (a.out/'input_checksums.json').write_text(json.dumps(checksums, indent=2)+'\n')
    report = {'status':'complete', 'python':platform.python_version(), 'command':sys.argv,
              'samples':samples, 'human_features':len(all_features),
              'verified_precursor_rows':sum(row['verified_precursor_rows'] for row in qc),
              'policy':'Unique counts exclude reads matching multiple distinct mature names; original multiplicities retained.',
              'limitations':['Conditional on existing alignments and references; no new alignment.',
                             'Fractional counts are sensitivity outputs, not integer count-model input.',
                             'Union counts still share reads across distinct mature names.',
                             'Viral counts retain historical BK-only reference attribution and BK/JC ambiguity.',
                             'No PCR deduplication or absolute molecule counting.']}
    (a.out/'run_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    (a.out/'SUCCESS').write_text('All sample reconstructions and original matrices verified.\n')
    print('SUCCESS: '+str(a.out),flush=True)


if __name__ == '__main__':
    main()
