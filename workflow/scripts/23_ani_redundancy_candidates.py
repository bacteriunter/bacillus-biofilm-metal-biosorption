#!/usr/bin/env python3
"""Summarize within-species skani sparse ANI pairs; never remove assemblies.

ANI threshold is fixed at 99.9%; AF 90% is a proposed, configurable threshold.
Connected components are candidate similarity groups, not proof that every
pair in a component passes the thresholds (single-linkage transitivity).
"""
import argparse
import csv
import json
import re
from collections import defaultdict, Counter
from pathlib import Path


def read_tsv(path):
    with path.open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def accession(path):
    m = re.search(r'(GC[AF]_\d+\.\d+)', Path(path).name)
    if not m:
        raise ValueError(f'Cannot identify accession from FASTA path: {path}')
    return m.group(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', type=Path, default=Path('data/processed/genome_cohort_selected.tsv'))
    ap.add_argument('--manifest', type=Path, default=Path('results/ani/module22/species_ani_manifest.tsv'))
    ap.add_argument('--sparse-dir', type=Path, default=Path('results/ani/module22/sparse'))
    ap.add_argument('--outdir', type=Path, default=Path('results/ani/module23'))
    ap.add_argument('--min-ani', type=float, default=99.9)
    ap.add_argument('--min-af-both', type=float, default=90.0)
    args = ap.parse_args()

    cohort = read_tsv(args.cohort)
    if not cohort:
        raise SystemExit('Empty cohort')
    by_acc = {r['assembly_accession']: r for r in cohort}
    if len(by_acc) != len(cohort):
        raise SystemExit('Duplicate assembly accessions in cohort')
    species_manifest = read_tsv(args.manifest)
    species_set = {r['species'] for r in species_manifest}
    if species_set != {r['species'] for r in cohort}:
        raise SystemExit('Species manifest and selected cohort disagree')
    if sum(int(r['n_genomes']) for r in species_manifest) != len(cohort):
        raise SystemExit('Species manifest counts do not sum to selected cohort')

    adjacency = {a: set() for a in by_acc}
    pairs = []
    expected_pairs = 0
    observed_pairs = 0
    missing = []
    for row in species_manifest:
        n = int(row['n_genomes'])
        if n < 2:
            continue
        expected_pairs += n * (n - 1) // 2
        source = Path(row['raw_ani_output'])
        path = args.sparse_dir / (source.stem + '.sparse.tsv')
        if not path.is_file():
            missing.append(str(path))
            continue
        with path.open(newline='', encoding='utf-8') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            required = {'Ref_file', 'Query_file', 'ANI', 'Align_fraction_ref', 'Align_fraction_query'}
            if not required.issubset(reader.fieldnames or []):
                raise SystemExit(f'Invalid sparse columns: {path}')
            seen = set()
            for r in reader:
                a, b = accession(r['Ref_file']), accession(r['Query_file'])
                if a == b or a not in by_acc or b not in by_acc:
                    raise SystemExit(f'Unexpected accession in {path}: {a}, {b}')
                if by_acc[a]['species'] != row['species'] or by_acc[b]['species'] != row['species']:
                    raise SystemExit(f'Cross-species pair in {path}: {a}, {b}')
                key = tuple(sorted((a, b)))
                if key in seen:
                    raise SystemExit(f'Duplicate pair in {path}: {key}')
                seen.add(key)
                ani = float(r['ANI'])
                af_ref = float(r['Align_fraction_ref'])
                af_query = float(r['Align_fraction_query'])
                qualifies = ani >= args.min_ani and min(af_ref, af_query) >= args.min_af_both
                if qualifies:
                    adjacency[a].add(b)
                    adjacency[b].add(a)
                pairs.append(dict(species=row['species'], assembly_a=a, assembly_b=b,
                                  ani_percent=ani, af_ref_percent=af_ref,
                                  af_query_percent=af_query, meets_proposed_thresholds='yes' if qualifies else 'no'))
            observed_pairs += len(seen)
    if missing:
        raise SystemExit('Missing sparse files (no outputs written):\n' + '\n'.join(missing[:20]))

    # Sparse skani may omit low-AF pairs; omitted pairs cannot be treated as observed nonmatches.
    components = []
    visited = set()
    for a in sorted(by_acc):
        if a in visited:
            continue
        stack = [a]
        visited.add(a)
        comp = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for b in sorted(adjacency[cur]):
                if b not in visited:
                    visited.add(b)
                    stack.append(b)
        components.append(sorted(comp))
    groups = []
    members = []
    for i, comp in enumerate(components, 1):
        group_id = f'ANI23_{i:04d}'
        priority = [a for a in comp if by_acc[a].get('phenotype_candidate', '').lower() == 'yes']
        groups.append(dict(group_id=group_id, species=by_acc[comp[0]]['species'],
                           n_genomes=len(comp), phenotype_candidate_count=len(priority),
                           accessions=';'.join(comp), phenotype_candidate_accessions=';'.join(priority),
                           interpretation='candidate_single_linkage_group' if len(comp)>1 else 'singleton'))
        for a in comp:
            members.append(dict(group_id=group_id, species=by_acc[a]['species'],
                                assembly_accession=a, phenotype_candidate=by_acc[a].get('phenotype_candidate',''),
                                n_group_members=len(comp), action='retain_pending_review'))

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir/'ani_pair_observations.tsv', pairs,
              ['species','assembly_a','assembly_b','ani_percent','af_ref_percent','af_query_percent','meets_proposed_thresholds'])
    write_tsv(args.outdir/'ani_candidate_groups.tsv', groups,
              ['group_id','species','n_genomes','phenotype_candidate_count','accessions','phenotype_candidate_accessions','interpretation'])
    write_tsv(args.outdir/'ani_candidate_members.tsv', members,
              ['group_id','species','assembly_accession','phenotype_candidate','n_group_members','action'])
    summary = dict(genomes=len(cohort), species=len(species_set), expected_within_species_pairs=expected_pairs,
                   reported_sparse_pairs=observed_pairs, omitted_sparse_pairs=expected_pairs-observed_pairs,
                   threshold_qualifying_pairs=sum(p['meets_proposed_thresholds']=='yes' for p in pairs),
                   candidate_groups=len(groups), multi_genome_groups=sum(len(c)>1 for c in components),
                   genomes_in_multi_genome_groups=sum(len(c) for c in components if len(c)>1),
                   min_ani_percent=args.min_ani, proposed_min_af_both_percent=args.min_af_both,
                   exclusions=0, limitations=['90% AF is proposed, not yet approved.',
                   'Single-linkage components may contain pairs not meeting both thresholds.',
                   'skani sparse omits low-alignment pairs; absence is not an observed ANI value.',
                   'No genomes selected for exclusion.'])
    (args.outdir/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
