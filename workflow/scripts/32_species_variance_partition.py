#!/usr/bin/env python3
"""Module 32: descriptive binary-variance partition by species (no inference)."""
import argparse, csv, hashlib, json
from collections import defaultdict
from pathlib import Path

def read(path):
    with path.open(newline='') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not reader.fieldnames: raise ValueError(f'Empty table: {path}')
        return list(reader),reader.fieldnames

def write(path,fields,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',default='data/processed/functional_cohort_638.tsv')
    p.add_argument('--matrix',default='results/functional/module28/gene_presence.tsv')
    p.add_argument('--outdir',default='results/functional/module32')
    a=p.parse_args()
    manifest,mcols=read(Path(a.manifest));matrix,cols=read(Path(a.matrix))
    if not {'assembly_accession','species'}<=set(mcols): raise ValueError('Manifest requires assembly_accession and species')
    if 'genome' not in cols:raise ValueError('Matrix requires genome')
    labels=[c for c in cols if c!='genome']
    by_acc={}
    for r in manifest:
        acc=r['assembly_accession'].strip();species=r['species'].strip()
        if not acc or not species:raise ValueError('Empty accession/species')
        if acc in by_acc and by_acc[acc]!=species:raise ValueError(f'Conflicting species: {acc}')
        by_acc[acc]=species
    by_species=defaultdict(list);seen=set()
    for r in matrix:
        genome=r['genome'].strip()
        if genome in seen:raise ValueError(f'Duplicate genome: {genome}')
        seen.add(genome)
        if genome not in by_acc:raise ValueError(f'Missing species for: {genome}')
        by_species[by_acc[genome]].append(r)
        for label in labels:
            if r[label] not in ('0','1'):raise ValueError(f'Nonbinary value {genome} {label}: {r[label]}')
    n=len(matrix)
    if n==0:raise ValueError('Empty matrix')
    output=[];species_rows=[]
    for label in labels:
        total_positive=sum(int(r[label]) for r in matrix)
        overall=total_positive/n
        # Exact identity for binary data, weighted by sampled genome counts:
        # p(1-p) = sum_s n_s/N * p_s(1-p_s) + sum_s n_s/N * (p_s-p)^2.
        within=0.;between=0.;variable_species=0;assessable=0
        for species,records in sorted(by_species.items()):
            ns=len(records);positive=sum(int(r[label]) for r in records);ps=positive/ns
            within+=(ns/n)*ps*(1-ps)
            between+=(ns/n)*(ps-overall)**2
            if ns>=2:
                assessable+=1
                variable_species+=int(0<positive<ns)
            species_rows.append({'group':label.split('::',1)[0] if '::' in label else '',
                                 'gene':label.split('::',1)[1] if '::' in label else label,
                                 'label':label,'species':species,'genomes':ns,
                                 'genomes_with_hit':positive,'prevalence':round(ps,8),
                                 'within_species_variable':int(ns>=2 and 0<positive<ns),
                                 'within_species_assessable':int(ns>=2)})
        total=overall*(1-overall)
        if abs(total-(within+between))>1e-10:raise AssertionError(f'Variance identity failed: {label}')
        output.append({'label':label,'genomes':n,'genomes_with_hit':total_positive,
                       'prevalence':round(overall,8),'total_binary_variance':round(total,10),
                       'between_species_component':round(between,10),
                       'within_species_component':round(within,10),
                       'between_fraction_of_total':round(between/total,8) if total else '',
                       'within_fraction_of_total':round(within/total,8) if total else '',
                       'species_with_2plus_genomes':assessable,
                       'species_variable_within':variable_species})
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    write(out/'gene_variance_partition.tsv',list(output[0]),output)
    write(out/'species_gene_detection.tsv',list(species_rows[0]),species_rows)
    summary={'status':'PASS','genomes':n,'species':len(by_species),'labels':len(labels),
             'variable_labels':sum(x['total_binary_variance']>0 for x in output),
             'species_with_2plus_genomes':sum(len(x)>=2 for x in by_species.values()),
             'single_genome_species':sum(len(x)==1 for x in by_species.values()),
             'method':'Genome-count-weighted binary variance decomposition: p(1-p) = weighted within-species p_s(1-p_s) + weighted between-species (p_s-p)^2.',
             'caution':'Descriptive sample partition, not phylogenetic control, statistical significance, validated orthology, or phenotype prediction.',
             'inputs':{str(path):hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in [a.manifest,a.matrix]},
             'outputs':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(out.glob('*.tsv'))}}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('inputs','outputs')},indent=2))
if __name__=='__main__':main()
