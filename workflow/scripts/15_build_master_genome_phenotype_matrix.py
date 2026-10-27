#!/usr/bin/env python3
"""Integrate existing bmetal evidence without treating co-occurrence as validation.
Standard library only; run from any directory.
"""
import csv
import json
from collections import defaultdict, Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / 'data/raw/pubmed_records.tsv'
NUM = ROOT / 'data/interim/pubmed_numeric_endpoint_triage.tsv'
MENTIONS = ROOT / 'data/interim/bacillus_taxon_mentions.tsv'
BIOSAMPLES = ROOT / 'data/interim/ncbi_biosample_resolution.tsv'
OUT = ROOT / 'data/processed/genome_phenotype_master.tsv'
LINKS = ROOT / 'data/processed/genome_phenotype_publication_links.tsv'
REPORT = ROOT / 'results/reports/genome_phenotype_master_summary.json'


def read_tsv(path, needed):
    if not path.is_file():
        raise SystemExit(f'ERROR: falta archivo: {path}')
    with path.open('r', encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        missing = set(needed) - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f'ERROR: columnas ausentes en {path.name}: {sorted(missing)}')
        return list(reader)


def write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def tokens(value):
    return {x.strip() for x in (value or '').split(';') if x.strip()}


def main():
    publications = read_tsv(RAW, ['pmid'])
    numeric = read_tsv(NUM, ['pmid', 'value_text', 'unit_text', 'context_text', 'triage_class'])
    mentions = read_tsv(MENTIONS, ['pmid', 'species_mention', 'strain_candidate', 'mention_level'])
    genomes = read_tsv(BIOSAMPLES, ['pmids', 'species_mention', 'strain_candidate', 'assembly_accession', 'biosample_accession', 'resolution_status'])

    known_pmids = {r['pmid'].strip() for r in publications if r['pmid'].strip()}
    mentions_by_pmid = defaultdict(set)
    for r in mentions:
        pmid = r['pmid'].strip()
        if pmid and r['strain_candidate'].strip():
            mentions_by_pmid[pmid].add((r['species_mention'].strip(), r['strain_candidate'].strip()))

    # Genome candidates are only associated with the PMID(s) recorded in module 09.
    genome_by_pmid = defaultdict(set)
    for r in genomes:
        for pmid in tokens(r['pmids']):
            genome_by_pmid[pmid].add((r['species_mention'].strip(), r['strain_candidate'].strip(),
                                      r['assembly_accession'].strip(), r['biosample_accession'].strip(),
                                      r['resolution_status'].strip()))

    # One observation per extracted numerical mention; do not replicate measurements across strains.
    columns = ['observation_id','pmid','doi','source_field','sentence_index','value_as_extracted',
               'unit_as_extracted','metal_mentions','phenotype_domain','triage_class',
               'evidence_text','species_in_measurement','strain_in_measurement',
               'biosorbent_material','dry_biomass_basis_confirmed','capacity_mg_g_dry',
               'experimental_conditions','phenotype_validation_status','phenotype_evidence_reference',
               'genome_link_status','genome_accession_verified','notes']
    master = []
    classes = Counter()
    for i, r in enumerate(numeric, 1):
        pmid = r['pmid'].strip()
        if not pmid:
            continue
        classes[r['triage_class']] += 1
        master.append({
            'observation_id': r.get('numeric_mention_id') or f'num_{i:06d}',
            'pmid':pmid, 'doi':r.get('doi',''), 'source_field':r.get('source_field',''),
            'sentence_index':r.get('sentence_index',''), 'value_as_extracted':r['value_text'],
            'unit_as_extracted':r['unit_text'], 'metal_mentions':r.get('metals',''),
            'phenotype_domain':r.get('domains',''), 'triage_class':r['triage_class'],
            'evidence_text':r['context_text'], 'species_in_measurement':'',
            'strain_in_measurement':'', 'biosorbent_material':'',
            'dry_biomass_basis_confirmed':'', 'capacity_mg_g_dry':'',
            'experimental_conditions':'', 'phenotype_validation_status':'unverified',
            'phenotype_evidence_reference':'', 'genome_link_status':'unverified',
            'genome_accession_verified':'',
            'notes':'Automated mention only; PMID-level co-occurrence cannot assign phenotype to strain.'
        })

    # Separate publication-level linkage table; no artificial Cartesian product with numeric values.
    link_columns = ['pmid','species_mention','strain_candidate','assembly_accession',
                    'biosample_accession','biosample_resolution_status','association_level',
                    'strain_phenotype_verified','genome_phenotype_verified']
    link_rows = []
    all_pmids = sorted(set(mentions_by_pmid) | set(genome_by_pmid))
    for pmid in all_pmids:
        genome_keys = set()
        for species, strain, assembly, biosample, status in sorted(genome_by_pmid[pmid]):
            genome_keys.add((species.casefold(), strain.casefold()))
            link_rows.append({'pmid':pmid,'species_mention':species,'strain_candidate':strain,
                              'assembly_accession':assembly,'biosample_accession':biosample,
                              'biosample_resolution_status':status,
                              'association_level':'publication_and_biosample_metadata_candidate',
                              'strain_phenotype_verified':'no', 'genome_phenotype_verified':'no'})
        for species, strain in sorted(mentions_by_pmid[pmid]):
            if (species.casefold(), strain.casefold()) not in genome_keys:
                link_rows.append({'pmid':pmid,'species_mention':species,'strain_candidate':strain,
                                  'assembly_accession':'','biosample_accession':'',
                                  'biosample_resolution_status':'',
                                  'association_level':'publication_strain_mention_only',
                                  'strain_phenotype_verified':'no', 'genome_phenotype_verified':'no'})

    write_tsv(OUT, columns, master)
    write_tsv(LINKS, link_columns, link_rows)
    report = {
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'literature_records':len(known_pmids),
        'numeric_mentions_in_master':len(master),
        'numeric_triage_classes':dict(sorted(classes.items())),
        'publications_with_numeric_mentions':len({r['pmid'] for r in master}),
        'publication_strain_genome_link_rows':len(link_rows),
        'publications_with_genome_candidates':len(genome_by_pmid),
        'verified_genome_phenotype_pairs':0,
        'master':str(OUT.relative_to(ROOT)),
        'publication_links':str(LINKS.relative_to(ROOT)),
        'methodological_note':'All original numeric classes preserved; no phenotype-to-strain assignment or genome-to-phenotype validation inferred. Missing genome in this table does not mean no genome exists.'
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(f'Matrix: {len(master)} numeric mentions; {len(link_rows)} publication-level strain/genome links.')
    print(f'Report: {REPORT.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
