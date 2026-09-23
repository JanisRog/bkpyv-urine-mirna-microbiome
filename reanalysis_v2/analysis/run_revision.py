#!/usr/bin/env python3
"""Rebuild current table analyses into a new output folder using explicit inputs."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

p = argparse.ArgumentParser(description=__doc__)
for name in ['profiles', 'mapping', 'metadata', 'recount', 'qc', 'out', 'rlib']:
    p.add_argument('--'+name, type=Path, required=True)
p.add_argument('--rscript', default='Rscript')
a = p.parse_args()
script = Path(__file__).resolve().parent
for name in ['profiles', 'mapping', 'metadata', 'recount', 'qc', 'rlib']:
    path = getattr(a, name).expanduser().resolve()
    if not path.exists():
        p.error(f'Missing {name}: {path}')
    setattr(a, name, path)
a.out = a.out.expanduser().resolve()
if a.out.exists():
    p.error('Output already exists. Use a new directory to preserve prior runs.')
a.out.mkdir(parents=True)

def run(argv):
    subprocess.run([str(x) for x in argv], check=True)

data = a.out / 'verified_data'
run([sys.executable, script/'rebuild_dna_tables.py', '--profiles', a.profiles,
     '--mapping', a.mapping, '--metadata', a.metadata, '--out', data])
run([a.rscript, script/'run_count_models.R', a.recount, a.metadata,
     a.out/'results/count_models', a.rlib])
run([sys.executable, script/'run_table_analyses.py', '--data', data, '--recount', a.recount,
     '--models', a.out/'results/count_models', '--qc', a.qc, '--out', a.out/'results/tables'])
run([sys.executable, script/'make_figures.py', '--root', a.out,
     '--recount', a.recount, '--metadata', a.metadata])
run([sys.executable, script/'validate_results.py', '--root', a.out, '--recount', a.recount])
versions = {'python': platform.python_version(), 'packages': {
    name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'scipy', 'matplotlib']}}
(a.out/'python_versions.json').write_text(json.dumps(versions, indent=2)+'\n')
(a.out/'SUCCESS').write_text('Table analyses and numerical validation completed; outstanding upstream DNA analyses remain separate.\n')
print('Completed:', a.out)
