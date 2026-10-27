#!/usr/bin/env python3
"""Compute skani ANI within each selected species, preserving raw outputs.

This module DOES NOT dereplicate, set aligned-fraction cutoffs, or discard genomes.
It produces per-species raw ANI matrices for a subsequent evidence-based decision.
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def fasta_ok(path):
    return path.is_file() and path.stat().st_size > 0 and path.open('rb').read(1) == b'>'


def slug_for(species):
    return re.sub(r'[^a-zA-Z0-9_-]+', '_', species).strip('_')[:70]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort', default='data/processed/genome_cohort_selected.tsv')
    ap.add_argument('--fastas', default='data/raw/genomes/cohort21')
    ap.add_argument('--output', default='results/ani/module22')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--max-species', type=int, default=0, help='For pilot only; 0 means all')
    ap.add_argument('--force', action='store_true', help='Recompute existing outputs')
    args = ap.parse_args()
    if not 1 <= args.threads <= 2:
        ap.error('Use 1 or 2 threads on this workstation')
    if args.max_species < 0:
        ap.error('--max-species must be >=0')
    if not shutil.which('skani'):
        sys.exit('ERROR: skani not found in PATH')
    with Path(args.cohort).open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or not {'assembly_accession','species','phenotype_candidate'}.issubset(rows[0]):
        sys.exit('ERROR: cohort TSV missing required columns')
    groups = defaultdict(list)
    seen = set()
    for row in rows:
        acc = row['assembly_accession'].strip()
        species = row['species'].strip()
        if not re.fullmatch(r'GCF_\d+\.\d+', acc) or not species or acc in seen:
            sys.exit(f'ERROR: invalid or repeated accession/species: {acc} / {species}')
        seen.add(acc)
        fasta = Path(args.fastas) / (acc + '.fna')
        if not fasta_ok(fasta):
            sys.exit(f'ERROR: missing/invalid FASTA {fasta}')
        groups[species].append((acc, fasta.resolve(), row['phenotype_candidate'].strip()))
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    species_sorted = sorted(groups)
    if args.max_species:
        species_sorted = species_sorted[:args.max_species]
    manifest = []
    for idx, species in enumerate(species_sorted, 1):
        items = sorted(groups[species])
        digest = hashlib.sha256(('\n'.join(f'{acc}\t{path}' for acc,path,_ in items)+'\n').encode()).hexdigest()[:12]
        stem = slug_for(species) + '_' + digest
        genome_list = output / (stem + '.list.txt')
        matrix = output / (stem + '.triangle.tsv')
        genome_list.write_text('\n'.join(str(path) for _,path,_ in items)+'\n', encoding='utf-8')
        status = 'singleton' if len(items) == 1 else 'existing'
        if len(items) > 1 and (args.force or not matrix.is_file() or matrix.stat().st_size == 0):
            print(f'[{idx}/{len(species_sorted)}] {species}: {len(items)} genomes', flush=True)
            tmp = output / (stem + '.triangle.partial')
            tmp.unlink(missing_ok=True)
            cmd = ['skani','triangle','-l',str(genome_list),'-o',str(tmp),'-t',str(args.threads)]
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode or not tmp.is_file() or tmp.stat().st_size == 0:
                tmp.unlink(missing_ok=True)
                print('ERROR executing:', ' '.join(cmd), file=sys.stderr)
                print((result.stderr or result.stdout)[-2000:], file=sys.stderr)
                sys.exit(f'ANI failed for {species}; completed matrices retained for resumption')
            tmp.replace(matrix)
            status = 'computed'
        manifest.append({'species':species,'n_genomes':len(items),'phenotype_candidate_genomes':sum(p=='yes' for _,_,p in items),'status':status,'genome_list':str(genome_list),'raw_ani_output':str(matrix) if len(items)>1 else '', 'accessions':';'.join(acc for acc,_,_ in items)})
    report = output / 'species_ani_manifest.tsv'
    with report.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=list(manifest[0]))
        w.writeheader(); w.writerows(manifest)
    print(json.dumps({'species_processed':len(manifest),'genomes_processed':sum(int(x['n_genomes']) for x in manifest),'species_with_comparisons':sum(x['n_genomes']>1 for x in manifest),'singleton_species':sum(x['n_genomes']==1 for x in manifest),'manifest':str(report),'note':'Raw skani results only. No genomes excluded or ANI clusters assigned.'},indent=2),flush=True)


if __name__ == '__main__':
    main()
