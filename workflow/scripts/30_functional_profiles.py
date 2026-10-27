#!/usr/bin/env python3
"""Module 30: descriptive functional-group profiles and pairwise co-detection.

Consumes validated module28 group presence and module29 genome richness. No new
sequence searches. Co-detection does not establish biological interaction.
"""
import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Empty TSV: {path}')
        return list(reader), reader.fieldnames


def write_tsv(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group-presence', default='results/functional/module28/group_presence.tsv')
    parser.add_argument('--genome-richness', default='results/functional/module29/genome_functional_richness.tsv')
    parser.add_argument('--outdir', default='results/functional/module30')
    args = parser.parse_args()
    group_path = Path(args.group_presence)
    richness_path = Path(args.genome_richness)
    groups_rows, group_columns = read_tsv(group_path)
    richness_rows, richness_columns = read_tsv(richness_path)
    groups = sorted(c for c in group_columns if c != 'genome')
    if not groups or 'genome' not in group_columns:
        raise ValueError('Group presence requires genome and group columns')
    if not {'genome', 'species'}.issubset(richness_columns):
        raise ValueError('Genome richness requires genome and species columns')
    if len(groups_rows) != len({r['genome'] for r in groups_rows}):
        raise ValueError('Duplicate genomes in group presence')
    if len(richness_rows) != len({r['genome'] for r in richness_rows}):
        raise ValueError('Duplicate genomes in genome richness')
    gr = {r['genome']: r for r in groups_rows}
    sp = {r['genome']: r['species'].strip() for r in richness_rows}
    if set(gr) != set(sp) or any(not x for x in sp.values()):
        raise ValueError('Genome sets or species labels do not match')
    states = {}
    for genome, row in gr.items():
        vals = {}
        for group in groups:
            value = row[group].strip()
            if value not in {'0', '1'}:
                raise ValueError(f'Nonbinary presence: {genome} {group}={value!r}')
            vals[group] = int(value)
        states[genome] = vals
    n = len(states)
    profile_counts = Counter()
    species_profile_counts = Counter()
    profile_rows = []
    for genome in sorted(states):
        active = [g for g in groups if states[genome][g]]
        profile = '|'.join(active) if active else 'none'
        profile_counts[profile] += 1
        species_profile_counts[(sp[genome], profile)] += 1
        profile_rows.append({'genome': genome, 'species': sp[genome], 'profile': profile,
                             'detected_groups': len(active), **states[genome]})
    profile_freq = [
        {'profile': p, 'genomes': count, 'prevalence': round(count / n, 6)}
        for p, count in sorted(profile_counts.items(), key=lambda x: (-x[1], x[0]))
    ]
    species_counts = Counter(sp.values())
    species_freq = [
        {'species': species, 'profile': profile, 'genomes_with_profile': count,
         'species_genomes': species_counts[species],
         'within_species_prevalence': round(count / species_counts[species], 6)}
        for (species, profile), count in sorted(species_profile_counts.items())
    ]
    pair_rows = []
    for a, b in itertools.combinations(groups, 2):
        n11 = sum(states[x][a] == 1 and states[x][b] == 1 for x in states)
        n10 = sum(states[x][a] == 1 and states[x][b] == 0 for x in states)
        n01 = sum(states[x][a] == 0 and states[x][b] == 1 for x in states)
        n00 = n - n11 - n10 - n01
        union = n11 + n10 + n01
        pair_rows.append({'group_a': a, 'group_b': b, 'both_present': n11,
                          'a_only': n10, 'b_only': n01, 'neither': n00,
                          'co_detection_prevalence': round(n11 / n, 6),
                          'jaccard': round(n11 / union, 6) if union else 0.0})
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / 'genome_group_profiles.tsv',
              ['genome', 'species', 'profile', 'detected_groups'] + groups, profile_rows)
    write_tsv(out / 'profile_frequencies.tsv', ['profile', 'genomes', 'prevalence'], profile_freq)
    write_tsv(out / 'species_profile_frequencies.tsv',
              ['species', 'profile', 'genomes_with_profile', 'species_genomes', 'within_species_prevalence'], species_freq)
    write_tsv(out / 'pairwise_group_codetection.tsv',
              ['group_a', 'group_b', 'both_present', 'a_only', 'b_only', 'neither',
               'co_detection_prevalence', 'jaccard'], pair_rows)
    summary = {
        'status': 'PASS', 'genomes': n, 'species': len(species_counts),
        'groups': groups, 'distinct_profiles': len(profile_counts),
        'pairwise_group_comparisons': len(pair_rows),
        'interpretation': ('Descriptive candidate-homolog co-detection only; not validated orthology, '
                           'gene function, biofilm formation, metal biosorption, or phenotype prediction.'),
        'inputs': {str(p): sha(p) for p in [group_path, richness_path]},
        'outputs': {p.name: sha(p) for p in sorted(out.glob('*.tsv'))}
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('inputs', 'outputs')}, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
