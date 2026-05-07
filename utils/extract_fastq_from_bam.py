#!/usr/bin/env python3
"""Extract paired FASTQs from SRA-downloaded BAMs via samtools collate | samtools fastq."""
import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',   required=True,      help='Metadata CSV with srr_accession (and optionally sample_id) columns')
    p.add_argument('--bamdir',  default='data/bam',  help='Directory containing {SRR}.bam files')
    p.add_argument('--outdir',  default='data/fastq', help='Output directory for FASTQ files')
    p.add_argument('--tmpdir',  default=None,         help='Temp prefix dir for samtools collate (default: outdir)')
    p.add_argument('--threads', type=int, default=6,  help='Threads for samtools (default: 6)')
    return p.parse_args()


def already_done(outdir: Path, srr: str) -> bool:
    return bool(list(outdir.glob(f'{srr}*.fastq')) or list(outdir.glob(f'{srr}*.fastq.gz')))


def fmt_elapsed(seconds: float) -> str:
    m, s = divmod(int(seconds), 60)
    return f'{m}m{s:02d}s' if m else f'{s}s'


def extract(bam: Path, srr: str, outdir: Path, tmpdir: Path, threads: int) -> bool:
    r1 = outdir / f'{srr}_1.fastq'
    r2 = outdir / f'{srr}_2.fastq'
    tmp_prefix = str(tmpdir / f'{srr}_collate_tmp')

    collate = subprocess.Popen(
        ['samtools', 'collate', '-O', '-u', '--no-PG',
         '-@', str(threads), str(bam), tmp_prefix],
        stdout=subprocess.PIPE,
    )
    fastq = subprocess.Popen(
        ['samtools', 'fastq',
         '-@', str(threads),
         '-1', str(r1),
         '-2', str(r2),
         '-0', '/dev/null',
         '-s', '/dev/null',
         '-n'],
        stdin=collate.stdout,
    )
    collate.stdout.close()
    fastq.communicate()
    collate.wait()

    if collate.returncode != 0 or fastq.returncode != 0:
        for f in (r1, r2):
            f.unlink(missing_ok=True)
        return False
    return True


def main():
    args = parse_args()
    bamdir  = Path(args.bamdir)
    outdir  = Path(args.outdir)
    tmpdir  = Path(args.tmpdir) if args.tmpdir else outdir
    outdir.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    with open(args.input) as fh:
        rows = [r for r in csv.DictReader(fh) if r.get('srr_accession')]

    total  = len(rows)
    failed = []
    t_start = time.time()

    for i, row in enumerate(rows, 1):
        srr = row['srr_accession']
        sid = row.get('sample_id', srr)
        prefix = f'[{i}/{total}]'
        bam = bamdir / f'{srr}.bam'

        if already_done(outdir, srr):
            print(f'{prefix} {sid} ({srr}) — skipping, FASTQ already exists', flush=True)
            continue

        if not bam.exists():
            print(f'{prefix} {sid} ({srr}) — skipping, BAM not found: {bam}', flush=True)
            continue

        qc = subprocess.run(['samtools', 'quickcheck', str(bam)], capture_output=True)
        if qc.returncode != 0:
            print(f'{prefix} {sid} ({srr}) — skipping, BAM incomplete or truncated', flush=True)
            continue

        print(f'\n{prefix} {sid} ({srr}) — extracting from {bam.name}', flush=True)
        t0 = time.time()

        ok = extract(bam, srr, outdir, tmpdir, args.threads)
        elapsed = time.time() - t0

        if ok:
            files = sorted(outdir.glob(f'{srr}*.fastq*'))
            sizes = ' + '.join(f'{f.stat().st_size / 1e9:.1f}G' for f in files)
            print(f'{prefix} {sid} ({srr}) — done in {fmt_elapsed(elapsed)} [{sizes}]', flush=True)
        else:
            print(f'{prefix} {sid} ({srr}) — FAILED', flush=True)
            failed.append((sid, srr))

    total_elapsed = fmt_elapsed(time.time() - t_start)
    print(f'\n{"=" * 60}')
    print(f'Finished {total - len(failed)}/{total} samples in {total_elapsed}')
    if failed:
        print(f'\nFailed ({len(failed)}):')
        for sid, srr in failed:
            print(f'  {sid}  {srr}')
        sys.exit(1)


if __name__ == '__main__':
    main()
