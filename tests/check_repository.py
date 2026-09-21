"""Check tracked files for common accidental data/path/output inclusions."""
import json
import re
import subprocess
from pathlib import Path

root=Path(__file__).resolve().parents[1]
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=root).decode().split('\0')
# Split literals so this checker does not flag its own source.
patterns=[r'/'+r'Users/[^/\s]+/',r'/'+r'home/[^/\s]+/',r'gh'+'[opusr]_[A-Za-z0-9]{20,}',r'github'+'_pat_[A-Za-z0-9_]{20,}']
allowed={'.py','.md','.txt','.ipynb','.yml','.yaml'}
fail=[]
for name in filter(None,tracked):
    p=root/name
    if p.suffix not in allowed and p.name!='.gitignore':fail.append(f'Unexpected tracked file type: {name}')
    text=subprocess.check_output(['git','show',':'+name],cwd=root).decode()
    if any(re.search(pattern,text) for pattern in patterns):fail.append(f'Private path or credential pattern: {name}')
    if p.suffix=='.ipynb':
        notebook=json.loads(text)
        for i,cell in enumerate(notebook['cells']):
            if cell.get('outputs') or cell.get('execution_count') is not None:
                fail.append(f'Notebook outputs/execution state: {name} cell {i}')
if fail:raise SystemExit('\n'.join(fail))
print(f'Checked {len(list(filter(None,tracked)))} tracked files; no flagged paths, credentials or notebook outputs.')
