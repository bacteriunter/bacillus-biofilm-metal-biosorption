#!/usr/bin/env python3
"""Module 24: auditable, non-destructive proposal of ANI representatives.

Only direct, observed pairs with ANI >= 99.9 and both AF >= 90 qualify.
Phenotype-candidate accessions are protected. No FASTA files are modified.
"""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ACCESSION = re.compile(r'GC[AF]_\d+\.\d+')
TRUE = {'yes', 'true', '1', 'y', 'si', 'sí'}

def load(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def save(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def number(row, fields, fallback):
    for key in fields:
        val = row.get(key, '').strip()
        if val:
            try:
                return float(val.replace('%', ''))
            except ValueError:
                pass
    return fallback

def rank(row):
    levels = {'complete genome': 4, 'chromosome': 3, 'scaffold': 2, 'contig': 1}
    level = next((row.get(k, '').strip().lower() for k in ('assembly_level', 'assemblyLevel', 'assembly_level_name') if row.get(k, '').strip()), '')
    completeness = number(row, ('completeness', 'checkm_completeness', 'checkm_completeness_percent'), -1)
    contamination = number(row, ('contamination', 'checkm_contamination', 'checkm_contamination_percent'), 1e9)
    return (levels.get(level, 0), completeness, -contamination)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', type=Path, default=Path('data/processed/genome_cohort_selected.tsv'))
    p.add_argument('--groups', type=Path, default=Path('results/ani/module23/ani_candidate_members.tsv'))
    p.add_argument('--pairs', type=Path, default=Path('results/ani/module23/ani_pair_observations.tsv'))
    p.add_argument('--phenotype-links', type=Path, default=Path('data/processed/phenotype_linked_subset.tsv'))
    p.add_argument('--outdir', type=Path, default=Path('results/ani/module24'))
    p.add_argument('--min-ani', type=float, default=99.9)
    p.add_argument('--min-af-both', type=float, default=90.0)
    args = p.parse_args()
    cohort = load(args.cohort)
    by_acc = {r['assembly_accession']: r for r in cohort}
    if len(by_acc) != len(cohort) or not by_acc:
        raise SystemExit('ERROR: duplicate accessions or empty cohort')
    members = load(args.groups)
    if len(members) != len(cohort) or {r['assembly_accession'] for r in members} != set(by_acc):
        raise SystemExit('ERROR: module23 members do not match selected cohort')
    group = {r['assembly_accession']: r['group_id'] for r in members}
    protected = {a for a, r in by_acc.items() if str(r.get('phenotype_candidate', '')).strip().lower() in TRUE}
    protected |= {r['assembly_accession'] for r in members if str(r.get('phenotype_candidate', '')).strip().lower() in TRUE}
    # Extract only exact assembly accessions from the phenotype-linked table.
    # Any ambiguous/non-assembly links remain candidates for manual review.
    phenotype_linked_exact = set()
    phenotype_link_rows = 0
    if args.phenotype_links.exists():
        for r in load(args.phenotype_links):
            phenotype_link_rows += 1
            for key, value in r.items():
                if value and any(word in key.lower() for word in ('accession', 'assembly', 'genome')):
                    phenotype_linked_exact.update(ACCESSION.findall(value))
        protected |= phenotype_linked_exact & set(by_acc)
    links_not_in_cohort = sorted(phenotype_linked_exact - set(by_acc))
    adj = {a: set() for a in by_acc}
    edge = {}
    observed = 0
    qualifying = 0
    for r in load(args.pairs):
        a, b = r['assembly_a'], r['assembly_b']
        if a not in by_acc or b not in by_acc or a == b:
            raise SystemExit(f'ERROR: invalid pair {a} {b}')
        if by_acc[a]['species'] != by_acc[b]['species']:
            raise SystemExit(f'ERROR: cross-species pair {a} {b}')
        key = tuple(sorted((a, b)))
        if key in edge:
            raise SystemExit(f'ERROR: duplicate pair {key}')
        ani = float(r['ani_percent'])
        af1, af2 = float(r['af_ref_percent']), float(r['af_query_percent'])
        good = ani >= args.min_ani and min(af1, af2) >= args.min_af_both
        edge[key] = (ani, af1, af2, good)
        observed += 1
        if good:
            qualifying += 1
            adj[a].add(b)
            adj[b].add(a)
    # Every proposed omission must have a direct qualifying edge to a retained genome.
    # Protect phenotype-linked accessions first; then favor better assembly quality.
    retained = set(protected)
    omissions = {}
    species_groups = defaultdict(list)
    for a, r in by_acc.items():
        species_groups[r['species']].append(a)
    for species in sorted(species_groups):
        order = sorted(species_groups[species], key=lambda a: (tuple(-v for v in rank(by_acc[a])), a))
        for a in order:
            if a in retained:
                continue
            possible = sorted((b for b in retained if b in adj[a]), key=lambda b: (tuple(-v for v in rank(by_acc[b])), b))
            if possible:
                rep = possible[0]
                ani, af1, af2, _ = edge[tuple(sorted((a, rep)))]
                omissions[a] = (rep, ani, af1, af2)
            else:
                retained.add(a)
    if retained & set(omissions) or retained | set(omissions) != set(by_acc) or not protected <= retained:
        raise SystemExit('ERROR: selection invariants failed')
    audit = []
    for a in sorted(by_acc, key=lambda x: (by_acc[x]['species'], x)):
        r = by_acc[a]
        candidate = a in protected
        om = omissions.get(a)
        audit.append({
            'assembly_accession': a,
            'species': r['species'],
            'module23_group_id': group[a],
            'phenotype_candidate_protected': 'yes' if candidate else 'no',
            'proposed_status': 'omit_redundant' if om else 'retain',
            'decision_reason': ('direct_ANI_AF_match_to_retained' if om else
                                'protected_phenotype_candidate' if candidate else 'no_qualifying_match_to_retained'),
            'representative_accession': om[0] if om else a,
            'direct_ani_percent': om[1] if om else '',
            'direct_af_ref_percent': om[2] if om else '',
            'direct_af_query_percent': om[3] if om else '',
            'assembly_level': r.get('assembly_level', ''),
            'completeness': r.get('completeness', r.get('checkm_completeness', '')),
            'contamination': r.get('contamination', r.get('checkm_contamination', '')),
        })
    columns = list(audit[0])
    args.outdir.mkdir(parents=True, exist_ok=True)
    save(args.outdir / 'representative_proposal_audit.tsv', audit, columns)
    save(args.outdir / 'proposed_retained.tsv', [r for r in audit if r['proposed_status'] == 'retain'], columns)
    save(args.outdir / 'proposed_omissions.tsv', [r for r in audit if r['proposed_status'] == 'omit_redundant'], columns)
    counts = Counter(r['species'] for r in audit if r['proposed_status'] == 'retain')
    summary = {
        'original_genomes': len(by_acc), 'species': len(species_groups),
        'observed_sparse_pairs': observed, 'qualifying_direct_pairs': qualifying,
        'proposed_retained_genomes': len(retained), 'proposed_redundant_omissions': len(omissions),
        'protected_phenotype_candidate_genomes': len(protected),
        'phenotype_link_table_rows': phenotype_link_rows,
        'phenotype_link_accessions_not_in_cohort': links_not_in_cohort,
        'species_with_retained_genomes': len(counts),
        'min_ani_percent': args.min_ani, 'min_af_both_percent': args.min_af_both,
        'files_deleted': 0, 'selection_applied_to_downstream': False,
        'method': 'greedy direct-pair representative proposal; phenotype accessions retained first; quality ranking where available',
        'limitations': [
            'Only observed sparse comparisons can support proposed omissions.',
            'Connected-component membership alone never justifies omission.',
            'Unverified phenotype links without exact assembly accessions require manual review.',
            'This is a proposal, not a finalized biological or analytical cohort.',
        ],
    }
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2, ensure_ascii=False))

if __name__ == '__main__':
    main()
