#!/usr/bin/env python3
"""Integrate literature evidence and provisional strain/genome matches by PMID.

A shared PMID establishes publication-level co-occurrence, NOT that a numeric
measurement belongs to a particular strain or genome.
"""
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EVIDENCE = ROOT / 'data/interim/pubmed_phenotype_evidence.tsv'
NUMERIC = ROOT / 'data/interim/pubmed_quantitative_mentions.tsv'
MENTIONS = ROOT / 'data/interim/bacillus_taxon_mentions.tsv'
GENOMES = ROOT / 'data/interim/ncbi_biosample_resolution.tsv'
OUT = ROOT / 'data/interim/pubmed_strain_phenotype_candidates.tsv'
SUMMARY = ROOT / 'results/reports/pubmed_strain_phenotype_integration_summary.json'


def read_tsv(path, required):
    if not path.is_file():
        raise SystemExit(f'Missing input: {path}')
    with path.open(encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f'Missing columns in {path}: {sorted(missing)}')
        return list(reader)


def pmid_list(value):
    return [s.strip() for s in (value or '').replace(',', ';').split(';') if s.strip()]


def main():
    ev = read_tsv(EVIDENCE, ('pmid','doi','source_field','sentence_index','evidence_text','domains','metals'))
    nums = read_tsv(NUMERIC, ('pmid','source_field','sentence_index','value_text','unit_text'))
    mentions = read_tsv(MENTIONS, ('pmid','species_mention','strain_candidate','mention_level'))
    genomes = read_tsv(GENOMES, ('pmids','species_mention','strain_candidate','assembly_accession','biosample_accession','resolution_status'))

    # Publication-level strain candidates: exclude genus/species-only mentions.
    strains_by_pmid = defaultdict(set)
    for r in mentions:
        if (r['strain_candidate'] or '').strip():
            strains_by_pmid[r['pmid'].strip()].add((r['species_mention'].strip(),r['strain_candidate'].strip()))

    # Preserve provenance of each possible genome, including multiple genomes per strain.
    genome_by_key = defaultdict(list)
    for r in genomes:
        for pmid in pmid_list(r['pmids']):
            key = (pmid, r['species_mention'].strip().casefold(), r['strain_candidate'].strip().casefold())
            genome_by_key[key].append(r)

    # Numeric quantities are joined to their own sentence, not arbitrarily to other passages.
    numbers_by_sentence = defaultdict(list)
    for r in nums:
        key = (r['pmid'].strip(),r['source_field'],r['sentence_index'])
        numbers_by_sentence[key].append(f"{r['value_text']} {r['unit_text']}")

    fields = ['pmid','doi','species_mention','strain_candidate','assembly_accession',
              'biosample_accession','genome_resolution_status','source_field','sentence_index',
              'domains','metals','numeric_mentions_in_sentence','evidence_text',
              'association_level','phenotype_measurement_status']
    OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    pub_ids = set()
    candidate_pub_ids = set()
    genome_pub_ids = set()
    with OUT.open('w',encoding='utf-8',newline='') as fh:
        writer = csv.DictWriter(fh,fieldnames=fields,delimiter='\t')
        writer.writeheader()
        for e in ev:
            pmid = e['pmid'].strip()
            pub_ids.add(pmid)
            candidates = sorted(strains_by_pmid.get(pmid, ()))
            if not candidates:
                counts['evidence_passages_without_strain_candidate'] += 1
                continue
            candidate_pub_ids.add(pmid)
            numbers = '; '.join(numbers_by_sentence.get((pmid,e['source_field'],e['sentence_index']), []))
            for species, strain in candidates:
                matches = genome_by_key.get((pmid,species.casefold(),strain.casefold()), [])
                if not matches:
                    matches = [None]
                else:
                    genome_pub_ids.add(pmid)
                for g in matches:
                    writer.writerow({
                        'pmid':pmid,'doi':e['doi'],'species_mention':species,'strain_candidate':strain,
                        'assembly_accession':g['assembly_accession'] if g else '',
                        'biosample_accession':g['biosample_accession'] if g else '',
                        'genome_resolution_status':g['resolution_status'] if g else 'no_biosample_candidate',
                        'source_field':e['source_field'],'sentence_index':e['sentence_index'],
                        'domains':e['domains'],'metals':e['metals'],
                        'numeric_mentions_in_sentence':numbers,'evidence_text':e['evidence_text'],
                        'association_level':'same_publication_only',
                        'phenotype_measurement_status':'not_verified_not_strain_assigned'})
                    counts['candidate_rows'] += 1
                    if numbers: counts['candidate_rows_with_numeric_mentions'] += 1
                    if g: counts['candidate_rows_with_biosample'] += 1
    report = {
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'input_evidence_rows':len(ev),'input_numeric_rows':len(nums),
        'input_strain_mention_rows':len(mentions),'input_biosample_resolution_rows':len(genomes),
        'publications_with_evidence':len(pub_ids),
        'publications_with_evidence_and_strain_candidate':len(candidate_pub_ids),
        'publications_with_evidence_strain_and_biosample_candidate':len(genome_pub_ids),
        **dict(sorted(counts.items())),
        'output':str(OUT.relative_to(ROOT)),
        'limitations':'All associations are publication-level only. No phenotype value has been attributed to a strain, no experimental endpoint validated, and no genome-phenotype link confirmed. Only genomes in module 09 are included.'
    }
    SUMMARY.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f"Integrated {len(ev)} passages; wrote {counts['candidate_rows']} candidate rows.")
    print(f'Report: {SUMMARY.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
