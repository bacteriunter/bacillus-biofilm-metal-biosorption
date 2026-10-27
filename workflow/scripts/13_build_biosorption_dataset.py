#!/usr/bin/env python3
"""Build an auditable biosorption evidence dataset; do not infer strain attribution."""
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRIAGE = ROOT / 'data/interim/pubmed_numeric_endpoint_triage.tsv'
LINK = ROOT / 'data/interim/pubmed_strain_phenotype_candidates.tsv'
OUT = ROOT / 'data/processed/biosorption_evidence_candidates.tsv'
GENOMES = ROOT / 'data/processed/biosorption_publication_genome_candidates.tsv'
REPORT = ROOT / 'results/reports/biosorption_dataset_build_summary.json'

def read(path, required):
    if not path.exists():
        raise SystemExit(f'Missing input: {path}')
    with path.open(encoding='utf-8-sig', newline='') as f:
        rd = csv.DictReader(f, delimiter='\t')
        missing = set(required) - set(rd.fieldnames or [])
        if missing:
            raise SystemExit(f'Missing columns in {path}: {sorted(missing)}')
        return list(rd)

def write(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)

def main():
    numeric = read(TRIAGE, ['pmid','doi','value_text','unit_text','context_text','triage_class','numeric_mention_id'])
    links = read(LINK, ['pmid','species_mention','strain_candidate','assembly_accession','biosample_accession'])
    candidates = [r for r in numeric if r['triage_class'] == 'candidate_metal_sorption_endpoint']
    # Preserve every candidate, including false positives, for explicit adjudication.
    evidence = []
    for r in candidates:
        evidence.append({
            'numeric_mention_id':r['numeric_mention_id'], 'pmid':r['pmid'], 'doi':r['doi'],
            'value_as_extracted':r['value_text'], 'unit_as_extracted':r['unit_text'],
            'metal_mentioned':r.get('metals',''), 'context_text':r['context_text'],
            'source_field':r.get('source_field',''), 'sentence_index':r.get('sentence_index',''),
            'organism_in_measurement':'', 'strain_in_measurement':'', 'biosorbent_material':'',
            'endpoint_type':'', 'reported_capacity_mg_g':'', 'experimental_conditions':'',
            'measurement_decision':'pending_fulltext_verification',
            'exclusion_reason':'', 'source_evidence_location':'',
            'genome_link_status':'not_evaluated',
            'notes':'Publication-level genomic candidates are in a separate table; do not assign them here.'
        })
    write(OUT, list(evidence[0]) if evidence else [
        'numeric_mention_id','pmid','doi','value_as_extracted','unit_as_extracted','metal_mentioned','context_text',
        'source_field','sentence_index','organism_in_measurement','strain_in_measurement','biosorbent_material',
        'endpoint_type','reported_capacity_mg_g','experimental_conditions','measurement_decision','exclusion_reason',
        'source_evidence_location','genome_link_status','notes'], evidence)
    pmids = {r['pmid'] for r in candidates}
    genome_rows = []
    seen = set()
    for r in links:
        if r['pmid'] not in pmids: continue
        key = (r['pmid'],r['species_mention'],r['strain_candidate'],r['assembly_accession'],r['biosample_accession'])
        if key in seen: continue
        seen.add(key)
        genome_rows.append({
            'pmid':r['pmid'],'species_mention':r['species_mention'],
            'strain_candidate':r['strain_candidate'],'assembly_accession':r['assembly_accession'],
            'biosample_accession':r['biosample_accession'],
            'genome_resolution_status':r.get('genome_resolution_status',''),
            'link_status':'same_publication_only_not_verified'
        })
    write(GENOMES, ['pmid','species_mention','strain_candidate','assembly_accession','biosample_accession','genome_resolution_status','link_status'], genome_rows)
    report = {
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'numeric_candidates':len(candidates),
        'unique_publications_with_numeric_candidates':len(pmids),
        'publications_with_genome_or_strain_context':len({r['pmid'] for r in genome_rows}),
        'publication_level_genome_context_rows':len(genome_rows),
        'evidence_table':str(OUT.relative_to(ROOT)),
        'genome_context_table':str(GENOMES.relative_to(ROOT)),
        'verified_phenotypes':0,
        'note':'No phenotype is validated or assigned to a strain. Full-text verification is required before modeling.'
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f"Prepared {len(candidates)} quantitative candidates across {len(pmids)} publications")
    print(f"Report: {REPORT.relative_to(ROOT)}")

if __name__ == '__main__':
    main()
