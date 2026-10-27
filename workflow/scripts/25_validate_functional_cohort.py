#!/usr/bin/env python3
"""Module 25: independently validate the approved provisional functional cohort.

Non-destructive. Does not rerun skani, alter original FASTA, or alter module24 outputs.
Run from the repository root: python workflow/scripts/25_validate_functional_cohort.py
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

ACC = re.compile(r'GC[AF]_\d+\.\d+')
TRUE = {'yes', 'true', '1', 'y', 'si', 'sí'}

def read_tsv(path):
    if not path.is_file():
        raise ValueError(f'Missing input: {path}')
    with path.open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Empty TSV: {path}')
        return list(reader), reader.fieldnames

def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', type=Path, default=Path('data/processed/genome_cohort_selected.tsv'))
    ap.add_argument('--audit', type=Path, default=Path('results/ani/module24/representative_proposal_audit.tsv'))
    ap.add_argument('--retained', type=Path, default=Path('results/ani/module24/proposed_retained.tsv'))
    ap.add_argument('--omissions', type=Path, default=Path('results/ani/module24/proposed_omissions.tsv'))
    ap.add_argument('--pairs', type=Path, default=Path('results/ani/module23/ani_pair_observations.tsv'))
    ap.add_argument('--phenotype-links', type=Path, default=Path('data/processed/phenotype_linked_subset.tsv'))
    ap.add_argument('--genomes-dir', type=Path, default=Path('data/raw/genomes/cohort21'))
    ap.add_argument('--outdir', type=Path, default=Path('results/validation/module25'))
    ap.add_argument('--min-ani', type=float, default=99.9)
    ap.add_argument('--min-af-both', type=float, default=90.0)
    ap.add_argument('--expected-original', type=int, default=877)
    ap.add_argument('--expected-retained', type=int, default=638)
    args = ap.parse_args()
    errors = []
    def check(ok, message):
        if not ok:
            errors.append(message)

    try:
        cohort, cohort_cols = read_tsv(args.cohort)
        audit, audit_cols = read_tsv(args.audit)
        kept, kept_cols = read_tsv(args.retained)
        omitted, omitted_cols = read_tsv(args.omissions)
        pairs, pair_cols = read_tsv(args.pairs)
        required = {
            'cohort': (cohort_cols, {'assembly_accession', 'species'}),
            'audit': (audit_cols, {'assembly_accession', 'species', 'proposed_status', 'representative_accession', 'phenotype_candidate_protected', 'direct_ani_percent', 'direct_af_ref_percent', 'direct_af_query_percent'}),
            'pairs': (pair_cols, {'assembly_a', 'assembly_b', 'ani_percent', 'af_ref_percent', 'af_query_percent'}),
        }
        for name, (columns, fields) in required.items():
            missing = fields - set(columns)
            if missing:
                raise ValueError(f'{name}: missing columns {sorted(missing)}')
    except (ValueError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1

    def index(rows, label):
        out = {}
        for r in rows:
            a = r.get('assembly_accession', '')
            if not a or a in out:
                check(False, f'{label}: empty or duplicate accession {a!r}')
            out[a] = r
        return out
    base = index(cohort, 'cohort')
    audit_by = index(audit, 'audit')
    kept_by = index(kept, 'retained')
    omitted_by = index(omitted, 'omissions')
    check(len(base) == args.expected_original, f'Expected {args.expected_original} original genomes; found {len(base)}')
    check(len(kept_by) == args.expected_retained, f'Expected {args.expected_retained} retained genomes; found {len(kept_by)}')
    check(set(base) == set(audit_by), 'Audit does not cover exactly the original cohort')
    check(set(kept_by).isdisjoint(omitted_by), 'Retained and omitted sets overlap')
    check(set(kept_by) | set(omitted_by) == set(base), 'Retained + omitted do not partition original cohort')
    check(len(kept) + len(omitted) == len(audit), 'Audit row counts inconsistent')

    observed = {}
    for p in pairs:
        a, b = p['assembly_a'], p['assembly_b']
        key = tuple(sorted((a, b)))
        if a == b or a not in base or b not in base or key in observed:
            check(False, f'Invalid or repeated observed pair {a} / {b}')
            continue
        try:
            vals = tuple(float(p[k]) for k in ('ani_percent', 'af_ref_percent', 'af_query_percent'))
        except (ValueError, TypeError):
            check(False, f'Invalid ANI/AF values for {a} / {b}')
            continue
        observed[key] = vals

    protected = {a for a, r in base.items() if str(r.get('phenotype_candidate', '')).strip().lower() in TRUE}
    protected |= {a for a, r in audit_by.items() if str(r.get('phenotype_candidate_protected', '')).strip().lower() in TRUE}
    exact_links = set()
    phenotype_rows = 0
    if args.phenotype_links.is_file():
        links, _ = read_tsv(args.phenotype_links)
        phenotype_rows = len(links)
        for row in links:
            for col, value in row.items():
                if value and any(word in col.lower() for word in ('accession', 'assembly', 'genome')):
                    exact_links.update(ACC.findall(value))
        protected |= exact_links & set(base)
    check(protected <= set(kept_by), f'Protected genomes omitted: {sorted(protected - set(kept_by))}')

    fasta_rows = []
    missing_fastas = []
    bad_fastas = []
    seen_paths = set()
    for accession in sorted(base):
        row = base[accession]
        arow = audit_by.get(accession, {})
        is_kept = accession in kept_by
        status = 'retain' if is_kept else 'omit_redundant'
        check(arow.get('proposed_status') == status, f'{accession}: audit status inconsistent')
        check(arow.get('species') == row['species'], f'{accession}: audit species mismatch')
        listed = kept_by.get(accession) if is_kept else omitted_by.get(accession)
        if listed:
            check(listed.get('proposed_status') == status, f'{accession}: subset status mismatch')
            check(listed.get('representative_accession') == arow.get('representative_accession'), f'{accession}: representative differs between audit and subset')
        rep = arow.get('representative_accession', '')
        if is_kept:
            check(rep == accession, f'{accession}: retained genome must represent itself')
        else:
            check(rep in kept_by, f'{accession}: representative {rep} is not retained')
            if rep in base:
                check(row['species'] == base[rep]['species'], f'{accession}: cross-species representative {rep}')
            edge = observed.get(tuple(sorted((accession, rep))))
            check(edge is not None, f'{accession}: no directly observed ANI pair to {rep}')
            if edge:
                ani, af1, af2 = edge
                check(ani >= args.min_ani and min(af1, af2) >= args.min_af_both,
                      f'{accession}: direct pair fails ANI/AF threshold ({ani}, {af1}, {af2})')
                for key, actual in zip(('direct_ani_percent', 'direct_af_ref_percent', 'direct_af_query_percent'), edge):
                    try:
                        check(abs(float(arow[key]) - actual) < 1e-6, f'{accession}: audit {key} differs from observed pair')
                    except (ValueError, TypeError, KeyError):
                        check(False, f'{accession}: missing/non-numeric audit {key}')
        path = args.genomes_dir / f'{accession}.fna'
        if not path.is_file():
            missing_fastas.append(accession)
        elif path.stat().st_size == 0:
            bad_fastas.append(accession)
        elif is_kept:
            try:
                with path.open('rb') as fh:
                    first = fh.read(1)
                if first != b'>':
                    bad_fastas.append(accession)
            except OSError:
                bad_fastas.append(accession)
        if is_kept:
            check(str(path) not in seen_paths, f'Duplicate FASTA path {path}')
            seen_paths.add(str(path))
            fasta_rows.append({
                'assembly_accession': accession,
                'species': row['species'],
                'fasta_path': str(path),
                'phenotype_candidate_protected': 'yes' if accession in protected else 'no',
                'assembly_level': row.get('assembly_level', ''),
                'completeness': row.get('completeness', row.get('checkm_completeness', '')),
                'contamination': row.get('contamination', row.get('checkm_contamination', '')),
                'module24_representative': rep,
            })
    check(not missing_fastas, f'Missing FASTA files ({len(missing_fastas)}): {missing_fastas[:12]}')
    check(not bad_fastas, f'Empty/invalid FASTA files ({len(bad_fastas)}): {bad_fastas[:12]}')
    check(len({r['species'] for r in fasta_rows}) == len({r['species'] for r in cohort}), 'Species lost from retained cohort')

    summary = {
        'validation_status': 'PASS' if not errors else 'FAIL',
        'original_genomes': len(base),
        'retained_genomes': len(kept_by),
        'proposed_omissions': len(omitted_by),
        'species_original': len({r['species'] for r in cohort}),
        'species_retained': len({r['species'] for r in fasta_rows}),
        'observed_direct_pairs': len(observed),
        'protected_phenotype_candidate_genomes': len(protected),
        'phenotype_link_table_rows': phenotype_rows,
        'phenotype_link_accessions_not_in_cohort': sorted(exact_links - set(base)),
        'missing_fasta_count': len(missing_fastas),
        'invalid_fasta_count': len(bad_fastas),
        'min_ani_percent': args.min_ani,
        'min_af_both_percent': args.min_af_both,
        'files_deleted': 0,
        'errors': errors[:100],
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / 'validation_summary.json'
    if errors:
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        print('FAILED: no downstream input manifest was generated', file=sys.stderr)
        return 1

    manifest = Path('data/processed/functional_cohort_638.tsv')
    write_tsv(manifest, fasta_rows, list(fasta_rows[0]))
    fasta_list = Path('data/processed/functional_cohort_638_fastas.txt')
    fasta_list.parent.mkdir(parents=True, exist_ok=True)
    fasta_list.write_text(''.join(r['fasta_path'] + '\n' for r in fasta_rows), encoding='utf-8')
    summary['functional_cohort_manifest'] = str(manifest)
    summary['functional_cohort_manifest_sha256'] = sha256(manifest)
    summary['fasta_list'] = str(fasta_list)
    summary['fasta_list_sha256'] = sha256(fasta_list)
    summary['input_sha256'] = {str(p): sha256(p) for p in (args.cohort, args.audit, args.retained, args.omissions, args.pairs)}
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0

if __name__ == '__main__':
    sys.exit(main())
