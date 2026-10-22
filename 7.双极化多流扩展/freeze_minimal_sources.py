"""Freeze this minimal validation's sources and inputs without altering prior manifests."""
from pathlib import Path
import os
import json
import hashlib

root = Path(__file__).resolve().parent
base = root.parent / '6.位置反馈方案'
files = [root / name for name in ['run_minimal_all_methods.m', 'verify_minimal_all_methods.m',
         'analyze_minimal_all_methods.py', 'build_minimal_reports.py', 'freeze_minimal_sources.py']]
files += list((root / '+min9').glob('*.m')) + list((root / '+dmimo').glob('*.m'))
for pattern in ['+upa/*.m', '+ula/*.m', '+posfb/*.m', 'lib/*.m', 'shared/*']:
    files += list(base.glob(pattern))
entries = [dict(path=os.path.relpath(p, root).replace('\\', '/'), bytes=p.stat().st_size,
                sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(files) if p.is_file()]
out = root / 'results/minimal_all_methods_20261022/source_and_input_SHA256.json'
out.write_text(json.dumps(dict(version='ALL9_MINIMAL_V1', files=entries), ensure_ascii=False, indent=2), encoding='utf-8')
print('Source/input file hashes:', len(entries))
