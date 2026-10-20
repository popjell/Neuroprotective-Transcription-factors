#!/usr/bin/env python3
"""Resume author Trobisch raw-count downloads; never replace existing inputs."""
from pathlib import Path
import subprocess
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scRNAseq/Trobisch_2022/data'
BASE = 'https://ms-cross-regional.cells.ucsc.edu/ms-cross-regional'
FILES = [(f'{r}/seurat.rds', f'{r}_seurat.rds') for r in ('cer', 'sc', 'cor')]
FILES += [('integrated/UMAP.coords.tsv.gz', 'UMAP.coords.tsv.gz')]
DATA.mkdir(parents=True, exist_ok=True)
for source, name in FILES:
    target = DATA / name
    if target.exists() and target.stat().st_size:
        print(f'Already present: {name}', flush=True)
        continue
    partial = target.with_name(name + '.part')
    print(f'Downloading {source}', flush=True)
    subprocess.run(['curl', '--fail', '--location', '--retry', '5', '--retry-delay', '5',
                    '--continue-at', '-', '--output', str(partial), f'{BASE}/{source}'], check=True)
    if not partial.stat().st_size:
        raise RuntimeError(f'Empty download: {source}')
    partial.rename(target)
    print(f'Complete: {name} ({target.stat().st_size:,} bytes)', flush=True)
