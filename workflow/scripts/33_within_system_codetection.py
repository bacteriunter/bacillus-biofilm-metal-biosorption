#!/usr/bin/env python3
"""Module 33: descriptive within-group candidate-homolog co-detection.

No alignments, inference of operons, or phenotype prediction.
"""
import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        return list(reader), reader.fieldnames


def write_tsv(path, fields, records):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene-presence', default='results/functional/module28/gene_presence.tsv')
    parser.add_argument('--manifest', default='data/processed/functional_cohort_638.tsv')
    parser.add_argument('--outdir', default='results/functional/module33')
    args = parser.parse_args()
    presence_path, manifest_path = Path(args.gene_presence), Path(args.manifest)
    output = Path(args.outdir)
    output.mkdir(parents=True, exist_ok=True)

    matrix, columns = read_tsv(presence_path)
    manifest, manifest_cols = read_tsv(manifest_path)
    required = {'assembly_accession', 'species'}
    if not required.issubset(manifest_cols):
        raise ValueError('Manifest must have assembly_accession and species')
    if not columns or columns[0] != 'genome':
        raise ValueError('Matrix first column must be genome')
    labels = columns[1:]
    if not labels or any('::' not in label for label in labels):
        raise ValueError('Expected group::gene columns')
    if len(set(labels)) != len(labels):
        raise ValueError('Duplicate labels')
    species_by_genome = {}
    for row in manifest:
        genome, species = row['assembly_accession'].strip(), row['species'].strip()
        if genome in species_by_genome and species_by_genome[genome] != species:
            raise ValueError(f'Conflicting species for {genome}')
        species_by_genome[genome] = species
    genomes = [row['genome'] for row in matrix]
    if len(set(genomes)) != len(genomes):
        raise ValueError('Duplicate genomes in matrix')
    if any(genome not in species_by_genome for genome in genomes):
        raise ValueError('Matrix genomes missing from manifest')
    for row in matrix:
        for label in labels:
            if row[label] not in ('0', '1'):
                raise ValueError(f'Nonbinary value for {row["genome"]}: {label}')

    by_group = defaultdict(list)
    for label in labels:
        group, gene = label.split('::', 1)
        by_group[group].append((gene, label))

    n = len(matrix)
    group_rows = []
    genome_rows = []
    pair_rows = []
    species_rows = []
    for group in sorted(by_group):
        genes = sorted(by_group[group])
        group_labels = [label for _, label in genes]
        full = sum(all(int(row[label]) for label in group_labels) for row in matrix)
        any_hit = sum(any(int(row[label]) for label in group_labels) for row in matrix)
        group_rows.append({'group': group, 'evaluated_gene_labels': len(genes),
                           'genomes': n, 'genomes_any_detected': any_hit,
                           'genomes_all_detected': full,
                           'fraction_all_detected': round(full/n, 8) if n else ''})
        for row in matrix:
            detected = sum(int(row[label]) for label in group_labels)
            genome_rows.append({'genome': row['genome'], 'species': species_by_genome[row['genome']],
                                'group': group, 'detected_labels': detected,
                                'evaluated_labels': len(group_labels),
                                'all_evaluated_detected': int(detected == len(group_labels))})
        for (gene_a, label_a), (gene_b, label_b) in itertools.combinations(genes, 2):
            both = sum(row[label_a] == '1' and row[label_b] == '1' for row in matrix)
            only_a = sum(row[label_a] == '1' and row[label_b] == '0' for row in matrix)
            only_b = sum(row[label_a] == '0' and row[label_b] == '1' for row in matrix)
            neither = n - both - only_a - only_b
            union = both + only_a + only_b
            pair_rows.append({'group': group, 'gene_a': gene_a, 'gene_b': gene_b,
                              'both': both, 'a_only': only_a, 'b_only': only_b,
                              'neither': neither, 'jaccard': round(both/union, 8) if union else ''})
        per_species = defaultdict(list)
        for row in matrix:
            per_species[species_by_genome[row['genome']]].append(row)
        for species, records in sorted(per_species.items()):
            all_count = sum(all(int(row[label]) for label in group_labels) for row in records)
            species_rows.append({'species': species, 'group': group,
                                 'genomes': len(records), 'genomes_all_evaluated_detected': all_count,
                                 'fraction_all_evaluated_detected': round(all_count/len(records), 8)})

    outputs = {
        'group_completeness_descriptive.tsv': (['group', 'evaluated_gene_labels', 'genomes', 'genomes_any_detected', 'genomes_all_detected', 'fraction_all_detected'], group_rows),
        'genome_group_detection_counts.tsv': (['genome', 'species', 'group', 'detected_labels', 'evaluated_labels', 'all_evaluated_detected'], genome_rows),
        'within_group_pairwise_codetection.tsv': (['group', 'gene_a', 'gene_b', 'both', 'a_only', 'b_only', 'neither', 'jaccard'], pair_rows),
        'species_group_all_detected.tsv': (['species', 'group', 'genomes', 'genomes_all_evaluated_detected', 'fraction_all_evaluated_detected'], species_rows),
    }
    for name, (fields, records) in outputs.items():
        write_tsv(output/name, fields, records)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {
        'status': 'PASS', 'genomes': n, 'species': len({species_by_genome[g] for g in genomes}),
        'gene_group_labels': len(labels), 'groups': len(by_group),
        'within_group_pairs': len(pair_rows),
        'caution': 'All evaluated labels detected does not establish complete operon, orthology, expression, biofilm or metal biosorption; reference panel may be incomplete.',
        'inputs': {str(path): digest(path) for path in (presence_path, manifest_path)},
        'outputs': {name: digest(output/name) for name in outputs}
    }
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({key: value for key, value in summary.items() if key not in ('inputs', 'outputs')}, indent=2))

if __name__ == '__main__':
    main()
