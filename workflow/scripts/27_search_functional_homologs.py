#!/usr/bin/env python3
"""Directed translated homology search for the validated Bacillus cohort.
Requires MMseqs2. No gene calling. Outputs are candidate homology alignments.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

FIELDS = 'query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,evalue,bits,qcov,tcov'

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta-list', default='data/processed/functional_cohort_638_fastas.txt')
    p.add_argument('--references', default='data/reference/functional_panel/functional_reference_panel_merged.faa')
    p.add_argument('--outdir', default='results/functional/module27')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--sensitivity', type=float, default=5.7)
    p.add_argument('--evalue', default='1e-5')
    p.add_argument('--force', action='store_true')
    args = p.parse_args()
    if not shutil.which('mmseqs'):
        p.error('mmseqs not found in PATH')
    refs = Path(args.references).resolve()
    if not refs.is_file():
        p.error(f'Reference FASTA missing: {refs}')
    paths = [Path(x.strip()).resolve() for x in Path(args.fasta_list).read_text().splitlines() if x.strip()]
    if len(paths) != 638 or len(set(paths)) != 638:
        p.error(f'Expected 638 unique genome paths; found {len(paths)}')
    missing = [str(x) for x in paths if not x.is_file()]
    if missing:
        p.error(f'Missing FASTA files: {missing[:5]}')
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'alignments').mkdir(exist_ok=True)
    (outdir / 'logs').mkdir(exist_ok=True)
    (outdir / 'tmp').mkdir(exist_ok=True)
    version = subprocess.check_output(['mmseqs', 'version'], text=True).strip()
    config = {'mmseqs_version': version, 'reference_sha256': sha256(refs),
              'fasta_list_sha256': sha256(args.fasta_list), 'fields': FIELDS,
              'search_type': 2, 'threads': args.threads, 'sensitivity': args.sensitivity,
              'evalue': args.evalue, 'genomes': len(paths)}
    config_path = outdir / 'run_config.json'
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        p.error('Configuration changed. Use a different --outdir to prevent mixing runs.')
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    manifest = outdir / 'progress.tsv'
    failures = []
    start = time.monotonic()
    with manifest.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['accession', 'status', 'alignments', 'seconds', 'genome_fasta'])
        for index, genome in enumerate(paths, 1):
            accession = genome.name.removesuffix('.fna')
            result = outdir / 'alignments' / f'{accession}.tsv'
            done = outdir / 'alignments' / f'{accession}.done.json'
            log = outdir / 'logs' / f'{accession}.log'
            tmp = outdir / 'tmp' / accession
            if done.exists() and result.exists() and not args.force:
                try:
                    metadata = json.loads(done.read_text())
                    if metadata.get('result_sha256') == sha256(result) and metadata.get('genome_sha256') == sha256(genome):
                        writer.writerow([accession, 'SKIPPED', metadata['alignments'], 0, genome])
                        handle.flush()
                        print(f'[{index}/638] {accession}: already complete', flush=True)
                        continue
                except (ValueError, KeyError, OSError):
                    pass
            done.unlink(missing_ok=True)
            result.unlink(missing_ok=True)
            shutil.rmtree(tmp, ignore_errors=True)
            partial = result.with_suffix('.partial.tsv')
            partial.unlink(missing_ok=True)
            cmd = ['mmseqs', 'easy-search', str(refs), str(genome), str(partial), str(tmp),
                   '--search-type', '2', '--threads', str(args.threads),
                   '-s', str(args.sensitivity), '-e', args.evalue,
                   '--format-output', FIELDS]
            t0 = time.monotonic()
            with log.open('w') as lf:
                lf.write('COMMAND: ' + ' '.join(cmd) + '\n')
                lf.flush()
                proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            elapsed = round(time.monotonic() - t0, 2)
            shutil.rmtree(tmp, ignore_errors=True)
            if proc.returncode != 0 or not partial.exists():
                failures.append(accession)
                partial.unlink(missing_ok=True)
                writer.writerow([accession, 'FAILED', '', elapsed, genome])
                print(f'[{index}/638] {accession}: FAILED (see {log})', flush=True)
                handle.flush()
                continue
            os.replace(partial, result)
            n = sum(1 for _ in result.open())
            done.write_text(json.dumps({'accession': accession, 'alignments': n,
                                        'result_sha256': sha256(result),
                                        'genome_sha256': sha256(genome),
                                        'run_config_sha256': sha256(config_path)}, indent=2) + '\n')
            writer.writerow([accession, 'COMPLETE', n, elapsed, genome])
            handle.flush()
            print(f'[{index}/638] {accession}: {n} alignments ({elapsed}s)', flush=True)
    summary = {'genomes': len(paths), 'failed': failures,
               'elapsed_seconds': round(time.monotonic() - start, 2),
               'output': str(outdir), 'interpretation': 'candidate homologs, not confirmed orthologs'}
    (outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    if failures:
        sys.exit(1)

if __name__ == '__main__':
    main()
