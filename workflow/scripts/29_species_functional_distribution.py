#!/usr/bin/env python3
"""Module 29: species-level descriptive prevalence from module 28, no alignments."""
import argparse, csv, json, hashlib
from collections import defaultdict, Counter
from pathlib import Path

def read_tsv(p):
    with open(p, newline='') as f:
        r=csv.DictReader(f, delimiter='\t')
        return list(r), r.fieldnames

def write_tsv(p, fields, rows):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def find_col(cols, names, label):
    lookup={c.lower().strip():c for c in cols}
    for n in names:
        if n in lookup:return lookup[n]
    raise ValueError(f'Cannot identify {label}; available columns: {cols}')

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',default='data/processed/genome_cohort_selected.tsv')
    p.add_argument('--gene-presence',default='results/functional/module28/gene_presence.tsv')
    p.add_argument('--group-presence',default='results/functional/module28/group_presence.tsv')
    p.add_argument('--outdir',default='results/functional/module29')
    a=p.parse_args(); out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    manifest, mc=read_tsv(a.manifest)
    accession_col=find_col(mc,['accession','assembly_accession','genome','genome_accession','assembly'], 'genome accession')
    species_col=find_col(mc,['species','species_name','organism_name','organism','taxon_name'], 'species')
    by_acc={}
    for row in manifest:
        acc=row[accession_col].strip()
        species=row[species_col].strip()
        if not acc or not species:continue
        if acc in by_acc and by_acc[acc]!=species:raise ValueError(f'Conflicting species for {acc}')
        by_acc[acc]=species
    genes,gcols=read_tsv(a.gene_presence)
    groups,pcols=read_tsv(a.group_presence)
    if not genes or not groups:raise ValueError('Module 28 matrices empty')
    gene_labels=[x for x in gcols if x!='genome'];group_labels=[x for x in pcols if x!='genome']
    gmap={r['genome']:r for r in genes};pmap={r['genome']:r for r in groups}
    if len(gmap)!=len(genes) or len(pmap)!=len(groups):raise ValueError('Duplicate genome in matrices')
    if set(gmap)!=set(pmap):raise ValueError('Gene and group matrices have different genome sets')
    unknown=set(gmap)-set(by_acc)
    if unknown:raise ValueError(f'Manifest lacks species for {len(unknown)} genomes: {sorted(unknown)[:5]}')
    species_to_genomes=defaultdict(list)
    for genome in sorted(gmap):species_to_genomes[by_acc[genome]].append(genome)
    gene_rows=[];group_rows=[];genome_rows=[];species_rows=[]
    for genome in sorted(gmap):
        genome_rows.append({'genome':genome,'species':by_acc[genome],
                            'detected_gene_labels':sum(int(gmap[genome][x]) for x in gene_labels),
                            'detected_groups':sum(int(pmap[genome][x]) for x in group_labels)})
    for species, genomes in sorted(species_to_genomes.items()):
        n=len(genomes)
        for label in gene_labels:
            count=sum(int(gmap[x][label]) for x in genomes)
            group,gene=label.split('::',1) if '::' in label else ('unassigned',label)
            gene_rows.append({'species':species,'group':group,'gene':gene,'genomes':n,
                              'genomes_with_hit':count,'prevalence':round(count/n,6)})
        for group in group_labels:
            count=sum(int(pmap[x][group]) for x in genomes)
            group_rows.append({'species':species,'group':group,'genomes':n,
                               'genomes_with_hit':count,'prevalence':round(count/n,6)})
        species_rows.append({'species':species,'genomes':n,
                             'gene_labels_detected_in_species':sum(any(int(gmap[x][c]) for x in genomes) for c in gene_labels),
                             'groups_detected_in_species':sum(any(int(pmap[x][c]) for x in genomes) for c in group_labels)})
    write_tsv(out/'species_gene_prevalence.tsv',['species','group','gene','genomes','genomes_with_hit','prevalence'],gene_rows)
    write_tsv(out/'species_group_prevalence.tsv',['species','group','genomes','genomes_with_hit','prevalence'],group_rows)
    write_tsv(out/'genome_functional_richness.tsv',['genome','species','detected_gene_labels','detected_groups'],genome_rows)
    write_tsv(out/'species_summary.tsv',['species','genomes','gene_labels_detected_in_species','groups_detected_in_species'],species_rows)
    summary={'status':'PASS','genomes':len(gmap),'species':len(species_to_genomes),
             'gene_labels':len(gene_labels),'groups':len(group_labels),
             'interpretation':'Descriptive prevalence of candidate homolog hits under module 28 thresholds; not confirmed orthology or biosorption.',
             'inputs':{str(x):hashlib.sha256(Path(x).read_bytes()).hexdigest() for x in [a.manifest,a.gene_presence,a.group_presence]},
             'outputs':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(out.glob('*.tsv'))}}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('inputs','outputs')},indent=2))
if __name__=='__main__':main()
