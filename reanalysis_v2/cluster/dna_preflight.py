#!/usr/bin/env python3
"""Read-only cluster inventory. Does not classify reads or change source files."""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def file_info(path):
    item = {'path': str(path), 'exists': path.is_file()}
    if item['exists']:
        st = path.stat()
        item.update(bytes=st.st_size, mtime_utc=datetime.datetime.fromtimestamp(
            st.st_mtime, datetime.timezone.utc).isoformat())
        # Avoid reading huge indexes, containers, alignments or read files.
        if st.st_size <= 10_000_000 and path.suffix in {'.py', '.sh', '.sbatch', '.txt', '.md'}:
            item['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return item


def fastq_prefix(path, limit=1000):
    lengths = collections.Counter()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        for _ in range(limit):
            header = handle.readline()
            if not header:
                break
            seq, plus, qual = (handle.readline().rstrip('\r\n') for _ in range(3))
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Malformed FASTQ prefix')
            lengths[len(seq)] += 1
    return {'records_examined': sum(lengths.values()), 'length_distribution': dict(lengths),
            'scope': 'First records only; not total yield or whole-file validation'}


def inventory(root, depth=3):
    rows = []
    if not root.is_dir():
        return [{'path': str(root), 'directory_missing': True}]
    for directory, dirs, files in os.walk(root, followlinks=False):
        if len(Path(directory).relative_to(root).parts) >= depth:
            dirs[:] = []
        for name in sorted(files):
            path = Path(directory) / name
            # Metadata only for database files; never parse core dumps.
            if name == 'core' or name.startswith('core.'):
                continue
            if (path.suffix in {'.fmi', '.bwt', '.sa', '.ann', '.amb', '.pac', '.bt2', '.bt2l', '.mmi', '.sif', '.img'}
                    or name in {'nodes.dmp', 'names.dmp'}
                    or path.suffix in {'.log', '.md'} or name.lower().startswith('readme')):
                rows.append(file_info(path))
    return rows


def command_record(argv):
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=45, check=False)
        return {'argv': argv, 'returncode': r.returncode,
                'output': (r.stdout + r.stderr)[:20000]}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {'argv': argv, 'error': str(exc)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, required=True)
    p.add_argument('--bioinfo', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not a.project.is_dir():
        p.error('Project directory does not exist')
    a.out.mkdir(parents=True, exist_ok=False)
    report = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'python': sys.version, 'project': str(a.project), 'samples': [], 'commands': []}
    for sample in sorted((a.project / 'results/metaphlan').glob('SRR*')):
        if not sample.is_dir():
            continue
        item = {'accession': sample.name, 'files': []}
        for suffix in ['_UC1.fastq', '_UC2.fastq', '_UC1.fastq.gz', '_UC2.fastq.gz',
                       '_con_sorted.bam', '_metaphlan4.bowtie2.bz2',
                       '_metaphlan4_profile_output.txt', '_con_alnstats.txt']:
            f = sample / (sample.name + suffix)
            item['files'].append(file_info(f))
            if f.is_file() and '.fastq' in suffix:
                try:
                    item['files'][-1]['prefix_check'] = fastq_prefix(f)
                except (OSError, ValueError, EOFError) as exc:
                    item['files'][-1]['prefix_error'] = str(exc)
        report['samples'].append(item)
    report['database_inventory'] = inventory(a.project / 'kaiju', 4)
    report['host_inventory'] = inventory(a.bioinfo / 'databases/Human_genome_Hg38', 1)
    report['metaphlan_database_inventory'] = inventory(a.project / 'metaphlan_db_local', 1)
    report['container_inventory'] = inventory(a.bioinfo / 'singularity', 1)
    report['script_inventory'] = [file_info(f) for f in sorted((a.project / '02_scripts').glob('*'))
                                  if f.is_file() and f.suffix in {'.sh', '.sbatch', '.py'}]
    engine = shutil.which('singularity') or shutil.which('apptainer')
    if engine:
        containers = a.bioinfo / 'singularity'
        for filename, toolargs in [
            ('metaphlan4.sif', ['metaphlan', '--version']),
            ('metaphlan4.sif', ['bowtie2', '--version']),
            ('samtools_1.15.1--h1170115_0.sif', ['samtools', '--version']),
            ('quay.io-biocontainers-bwa-0.7.17--h5bf99c6_8.img', ['bwa']),
        ]:
            if (containers / filename).is_file():
                report['commands'].append(command_record([engine, 'exec', str(containers / filename), *toolargs]))
    candidates = [Path.home() / 'data/conda/envs/kaiju/bin/kaiju']
    if shutil.which('kaiju'):
        candidates.append(Path(shutil.which('kaiju')))
    for binary in dict.fromkeys(candidates):
        if binary.is_file():
            report['commands'].append(command_record([str(binary), '-h']))
    report['limitations'] = [
        'Current executable versions do not independently establish versions used historically.',
        'Database filenames and modification times do not establish reference release or taxonomic coverage.',
        'This inventory does not measure residual human DNA or classify viruses, fungi or parasites.',
        'FASTQ examination is prefix-only. No sequence strings or read identifiers are written to this report.',
    ]
    (a.out / 'dna_preflight.json').write_text(json.dumps(report, indent=2) + '\n')
    (a.out / 'SUCCESS').write_text('Inventory completed. Review report for missing files and tool errors.\n')
    print('Report:', a.out / 'dna_preflight.json')


if __name__ == '__main__':
    main()
