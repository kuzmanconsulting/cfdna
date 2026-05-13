#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description='Download cfDNA BAMs via sam-dump | samtools')
    p.add_argument('--input',   required=True,     help='Metadata CSV with srr_accession column')
    p.add_argument('--outdir',  default='data/bam', help='Output directory for BAM files')
    p.add_argument('--threads', type=int, default=6, help='Threads for samtools view (default: 6)')
    return p.parse_args()


def already_done(outdir: Path, srr: str) -> bool:
    return (outdir / f'{srr}.bam').exists()


def fmt_elapsed(seconds: float) -> str:
    m, s = divmod(int(seconds), 60)
    return f'{m}m{s:02d}s' if m else f'{s}s'


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.input) as fh:
        rows = [r for r in csv.DictReader(fh) if r.get('srr_accession')]

    total = len(rows)
    failed = []
    t_start = time.time()

    for i, row in enumerate(rows, 1):
        sid = row['sample_id']
        srr = row['srr_accession']
        prefix = f'[{i}/{total}]'
        out_bam = outdir / f'{srr}.bam'

        if already_done(outdir, srr):
            print(f'{prefix} {sid} ({srr}) — skipping, BAM already exists', flush=True)
            continue

        print(f'\n{prefix} {sid} ({srr}) — starting', flush=True)
        t0 = time.time()

        sam_dump = subprocess.Popen(
            ['sam-dump', '--header', srr],
            stdout=subprocess.PIPE,
        )
        samtools = subprocess.Popen(
            ['samtools', 'view', '-bS', '-@', str(args.threads), '-o', str(out_bam)],
            stdin=sam_dump.stdout,
        )
        sam_dump.stdout.close()
        samtools.communicate()
        sam_dump.wait()
        elapsed = time.time() - t0

        if sam_dump.returncode == 0 and samtools.returncode == 0:
            size = out_bam.stat().st_size / 1e9
            print(f'{prefix} {sid} ({srr}) — done in {fmt_elapsed(elapsed)} [{size:.1f}G]', flush=True)
        else:
            print(f'{prefix} {sid} ({srr}) — FAILED (sam-dump={sam_dump.returncode} samtools={samtools.returncode})', flush=True)
            out_bam.unlink(missing_ok=True)
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
