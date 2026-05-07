#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description='Download cfDNA FASTQs via fasterq-dump')
    p.add_argument('--input',   required=True,      help='Metadata CSV with srr_accession column')
    p.add_argument('--outdir',  default='data/fastq',help='Output directory for FASTQ files')
    p.add_argument('--threads', type=int, default=6, help='Threads per fasterq-dump call (default: 6)')
    p.add_argument('--tmpdir',  default=None,        help='Temp dir for fasterq-dump (default: outdir)')
    return p.parse_args()


def already_done(outdir: Path, srr: str) -> bool:
    return bool(list(outdir.glob(f'{srr}*.fastq')) or list(outdir.glob(f'{srr}*.fastq.gz')))


def fmt_elapsed(seconds: float) -> str:
    m, s = divmod(int(seconds), 60)
    return f'{m}m{s:02d}s' if m else f'{s}s'


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(args.tmpdir) if args.tmpdir else outdir

    with open(args.input) as fh:
        rows = [r for r in csv.DictReader(fh) if r.get('srr_accession')]

    total = len(rows)
    failed = []
    t_start = time.time()

    for i, row in enumerate(rows, 1):
        sid = row['sample_id']
        srr = row['srr_accession']
        prefix = f'[{i}/{total}]'

        if already_done(outdir, srr):
            print(f'{prefix} {sid} ({srr}) — skipping, FASTQ already exists', flush=True)
            continue

        print(f'\n{prefix} {sid} ({srr}) — starting', flush=True)
        t0 = time.time()

        cmd = [
            'fasterq-dump', srr,
            '--split-files',
            '--outdir', str(outdir),
            '--temp',   str(tmpdir),
            '--threads', str(args.threads),
        ]

        result = subprocess.run(cmd)
        elapsed = time.time() - t0

        if result.returncode == 0:
            files = sorted(outdir.glob(f'{srr}*.fastq*'))
            sizes = ' + '.join(f'{f.stat().st_size/1e9:.1f}G' for f in files)
            print(f'{prefix} {sid} ({srr}) — done in {fmt_elapsed(elapsed)} [{sizes}]', flush=True)
        else:
            print(f'{prefix} {sid} ({srr}) — FAILED (exit {result.returncode})', flush=True)
            failed.append((sid, srr))

    total_elapsed = fmt_elapsed(time.time() - t_start)
    print(f'\n{"="*60}')
    print(f'Finished {total - len(failed)}/{total} samples in {total_elapsed}')
    if failed:
        print(f'\nFailed ({len(failed)}):')
        for sid, srr in failed:
            print(f'  {sid}  {srr}')
        sys.exit(1)


if __name__ == '__main__':
    main()
