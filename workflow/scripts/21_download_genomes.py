#!/usr/bin/env python3
"""Download selected RefSeq genomic FASTAs in resumable, verified NCBI Datasets batches.

Only downloads genomes from the module 20 selected TSV. Does not run ANI.
Python >=3.10, NCBI datasets CLI. No third-party Python dependencies.
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ACCESSION = re.compile(r'^GCF_\d+\.\d+$')


def read_rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or 'assembly_accession' not in rows[0]:
        raise ValueError('Selected TSV empty or missing assembly_accession')
    seen = set()
    for row in rows:
        acc = row['assembly_accession'].strip()
        if not ACCESSION.fullmatch(acc) or acc in seen:
            raise ValueError(f'Invalid or repeated assembly accession: {acc}')
        seen.add(acc)
    return rows


def valid_fasta(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    with path.open('rb') as handle:
        return handle.read(1) == b'>'


def extract_batch(archive, batch, output):
    """Return accession->reason for each missing or invalid FASTA in archive."""
    failures = {}
    with zipfile.ZipFile(archive) as zf:
        members = {}
        for info in zf.infolist():
            if info.is_dir():
                continue
            m = re.search(r'(?:^|/)ncbi_dataset/data/(GCF_\d+\.\d+)/[^/]+\.fna$', info.filename)
            if m:
                members.setdefault(m.group(1), []).append(info)
        for acc in batch:
            candidates = members.get(acc, [])
            if len(candidates) != 1:
                failures[acc] = f'expected_one_fna_found_{len(candidates)}'
                continue
            destination = output / f'{acc}.fna'
            tmp = output / f'{acc}.fna.partial'
            try:
                with zf.open(candidates[0]) as source, tmp.open('wb') as target:
                    shutil.copyfileobj(source, target)
                if not valid_fasta(tmp):
                    raise ValueError('empty_or_invalid_fasta')
                tmp.replace(destination)
            except Exception as exc:
                tmp.unlink(missing_ok=True)
                failures[acc] = str(exc)
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', default='data/processed/genome_cohort_selected.tsv')
    ap.add_argument('--output', default='data/raw/genomes/cohort21')
    ap.add_argument('--batch-size', type=int, default=25)
    ap.add_argument('--max-batches', type=int, default=0, help='0 means all; use 1 for pilot')
    ap.add_argument('--retries', type=int, default=2)
    args = ap.parse_args()
    if args.batch_size < 1 or args.max_batches < 0 or args.retries < 0:
        ap.error('batch-size >=1, max-batches >=0, retries >=0 required')
    if not shutil.which('datasets'):
        sys.exit('ERROR: datasets CLI not found')
    rows = read_rows(Path(args.manifest))
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    all_accessions = [r['assembly_accession'].strip() for r in rows]
    missing = [acc for acc in all_accessions if not valid_fasta(output / f'{acc}.fna')]
    batches = [missing[i:i + args.batch_size] for i in range(0, len(missing), args.batch_size)]
    if args.max_batches:
        batches = batches[:args.max_batches]
    print(f'Selected: {len(all_accessions)} | already present: {len(all_accessions)-len(missing)} | missing: {len(missing)} | batches this run: {len(batches)}', flush=True)
    errors = {}
    for idx, batch in enumerate(batches, 1):
        print(f'Batch {idx}/{len(batches)}: {len(batch)} accessions', flush=True)
        batch_errors = {acc: 'download_not_completed' for acc in batch}
        for attempt in range(1, args.retries + 2):
            pending = [acc for acc in batch if not valid_fasta(output / f'{acc}.fna')]
            if not pending:
                break
            with tempfile.TemporaryDirectory(prefix='ncbi_batch_') as tmpdir:
                tmp = Path(tmpdir)
                accession_file = tmp / 'accessions.txt'
                accession_file.write_text('\n'.join(pending) + '\n', encoding='utf-8')
                archive = tmp / 'genomes.zip'
                command = ['datasets', 'download', 'genome', 'accession', '--inputfile', str(accession_file), '--include', 'genome', '--filename', str(archive)]
                result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
                if result.returncode or not archive.exists():
                    msg = (result.stderr or 'no zip produced').strip()[-600:]
                    batch_errors.update({acc: f'datasets_failed: {msg}' for acc in pending})
                    print(f'  attempt {attempt} failed: {msg}', file=sys.stderr, flush=True)
                    continue
                try:
                    failures = extract_batch(archive, pending, output)
                    batch_errors = failures
                except (zipfile.BadZipFile, OSError) as exc:
                    batch_errors = {acc: f'zip_error: {exc}' for acc in pending}
                    print(f'  attempt {attempt} extraction failed: {exc}', file=sys.stderr, flush=True)
        for acc in batch:
            if not valid_fasta(output / f'{acc}.fna'):
                errors[acc] = batch_errors.get(acc, 'missing_after_download')
    status = []
    for row in rows:
        acc = row['assembly_accession'].strip()
        path = output / f'{acc}.fna'
        status.append({'assembly_accession': acc, 'species': row.get('species', ''), 'phenotype_candidate': row.get('phenotype_candidate', ''), 'fasta_path': str(path) if valid_fasta(path) else '', 'status': 'downloaded' if valid_fasta(path) else 'missing', 'error': errors.get(acc, '')})
    report = output / 'download_manifest.tsv'
    with report.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(status[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(status)
    complete = sum(r['status'] == 'downloaded' for r in status)
    print(json.dumps({'selected': len(status), 'downloaded': complete, 'remaining': len(status) - complete, 'manifest': str(report)}, indent=2), flush=True)
    if errors:
        print(f'Failed accessions this run: {len(errors)}; rerun to retry', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
