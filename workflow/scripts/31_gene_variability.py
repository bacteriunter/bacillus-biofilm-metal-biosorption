#!/usr/bin/env python3
"""Module 31: descriptive gene-label prevalence and between-species heterogeneity.
Uses module28 binary gene presence and validated 638-genome manifest.
No alignment, phenotype inference, or statistical hypothesis testing.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with path.open(newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        if not r.fieldnames:
            raise ValueError(f'Empty TSV: {path}')
        return list(r), r.fieldnames


def write_tsv(path, cols, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default='data/processed/functional_cohort_638.tsv')
    ap.add_argument('--gene-presence', default='results/functional/module28/gene_presence.tsv')
    ap.add_argument('--outdir', default='results/functional/module31')
    a = ap.parse_args()
    manifest_path = Path(a.manifest)
    matrix_path = Path(a.gene_presence)
    manifest, mcols = read_tsv(manifest_path)
    matrix, cols = read_tsv(matrix_path)
    if not {'assembly_accession', 'species'} <= set(mcols):
        raise ValueError('Manifest must contain assembly_accession and species')
    if 'genome' not in cols:
        raise ValueError('Gene matrix must contain genome')
    labels = [c for c in cols if c != 'genome']
    if len(labels) != len(set(labels)) or any('::' not in c for c in labels):
        raise ValueError('Duplicate or malformed group::gene labels')
    species_of = {}
    for r in manifest:
        acc = r['assembly_accession'].strip()  # Preserve accession version, e.g. GCF_... .1
        species = r['species'].strip()
        if not acc or not species:
            raise ValueError('Blank accession/species in manifest')
        if acc in species_of:
            raise ValueError(f'Duplicate accession in manifest: {acc}')
        species_of[acc] = species
    presence = {}
    for r in matrix:
        genome = r['genome'].strip()
        if genome in presence:
            raise ValueError(f'Duplicate genome in gene matrix: {genome}')
        try:
            presence[genome] = {c: int(r[c]) for c in labels}
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError(f'Invalid presence value for {genome}: {exc}') from exc
        if any(v not in (0, 1) for v in presence[genome].values()):
            raise ValueError(f'Nonbinary presence values for {genome}')
    if set(presence) != set(species_of):
        missing = sorted(set(presence) - set(species_of))[:5]
        extra = sorted(set(species_of) - set(presence))[:5]
        raise ValueError(f'Genome mismatch; missing species {missing}; extra manifest {extra}')

    species_genomes = defaultdict(list)
    for genome in sorted(presence):
        species_genomes[species_of[genome]].append(genome)
    n = len(presence)
    gene_rows = []
    species_rows = []
    species_label_rows = []
    for species, genomes in sorted(species_genomes.items()):
        ns = len(genomes)
        variable = 0
        for label in labels:
            count = sum(presence[g][label] for g in genomes)
            if 0 < count < ns:
                variable += 1
            group, gene = label.split('::', 1)
            species_label_rows.append({
                'species': species, 'genomes': ns, 'group': group, 'gene': gene,
                'genomes_with_hit': count, 'prevalence': round(count / ns, 6),
                'within_species_variable': int(0 < count < ns),
                'within_species_assessable': int(ns >= 2),
            })
        species_rows.append({
            'species': species, 'genomes': ns,
            'within_species_assessable': int(ns >= 2),
            'variable_gene_labels_within_species': variable if ns >= 2 else '',
        })

    for label in labels:
        group, gene = label.split('::', 1)
        count = sum(presence[g][label] for g in presence)
        present_species = 0
        absent_species = 0
        variable_species = 0
        for genomes in species_genomes.values():
            k = sum(presence[g][label] for g in genomes)
            present_species += int(k > 0)
            absent_species += int(k == 0)
            variable_species += int(0 < k < len(genomes))
        gene_rows.append({
            'group': group, 'gene': gene, 'genomes': n,
            'genomes_with_hit': count, 'genome_prevalence': round(count / n, 6),
            'status': 'universal' if count == n else ('undetected' if count == 0 else 'variable'),
            'species_total': len(species_genomes),
            'species_with_any_hit': present_species,
            'species_without_hit': absent_species,
            'species_within_species_variable': variable_species,
        })

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        'gene_global_prevalence.tsv': (['group','gene','genomes','genomes_with_hit','genome_prevalence','status','species_total','species_with_any_hit','species_without_hit','species_within_species_variable'], gene_rows),
        'species_gene_prevalence.tsv': (['species','genomes','group','gene','genomes_with_hit','prevalence','within_species_variable','within_species_assessable'], species_label_rows),
        'species_variability_summary.tsv': (['species','genomes','within_species_assessable','variable_gene_labels_within_species'], species_rows),
    }
    for filename, (headers, data) in paths.items():
        write_tsv(out / filename, headers, data)
    summary = {
        'status': 'PASS', 'genomes': n, 'species': len(species_genomes),
        'gene_group_labels': len(labels),
        'universal_labels': sum(r['status'] == 'universal' for r in gene_rows),
        'variable_labels': sum(r['status'] == 'variable' for r in gene_rows),
        'undetected_labels': sum(r['status'] == 'undetected' for r in gene_rows),
        'species_with_at_least_two_genomes': sum(len(gs) >= 2 for gs in species_genomes.values()),
        'single_genome_species': sum(len(gs) == 1 for gs in species_genomes.values()),
        'interpretation': 'Descriptive candidate homolog detection; not validated orthology, phenotype or biosorption. Species with one genome do not support within-species variability assessment.',
        'inputs': {str(p): sha(p) for p in (manifest_path, matrix_path)},
        'outputs': {name: sha(out / name) for name in paths},
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('inputs','outputs')}, indent=2))

if __name__ == '__main__':
    main()
