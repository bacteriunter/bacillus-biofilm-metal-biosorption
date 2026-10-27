#!/usr/bin/env python3
"""Conservative, auditable title/abstract triage of bmetal PubMed records.

This is prioritization, not a final eligibility assessment.
No records are discarded; manual review is mandatory for ambiguous cases.
"""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

BIOFILM = re.compile(r'\bbiofilms?\b|\bextracellular matrix\b|\bextracellular polymeric substances\b|\bEPS\b', re.I)
METAL_PROCESS = re.compile(r'\bbiosorp\w*\b|\bbioadsorp\w*\b|\bmetal (?:adsorption|binding|uptake)\b|\badsorp\w*\b', re.I)
METAL = re.compile(r'\b(?:metal|metals|cadmium|lead|chromium|copper|nickel|zinc|arsenic|mercury|Cd\(II\)|Pb\(II\)|Cr\(VI\))\b', re.I)
QUANT = re.compile(r'\b(?:quantif\w*|measur\w*|assay|assays|crystal violet|absorbance|optical density|mg/g|mg g-1|removal efficiency|adsorption capacity|percent|percentage|%)\b', re.I)
REVIEW_TITLE = re.compile(r'\b(?:systematic review|scoping review|narrative review|literature review|meta-analysis|bibliometric (?:analysis|review))\b|\b(?:a|an) review\b|\breview of\b', re.I)
FIELDS = ['pmid','doi','title','authors','year','journal','abstract','search_groups']
EXTRA = ['screen_category','screen_reason','biofilm_signal','biosorption_signal','quantitative_signal','screen_method']

def screen(rec):
    title = rec['title']
    combined = title + ' ' + rec['abstract']
    bio = bool(BIOFILM.search(combined))
    sorp = bool(METAL_PROCESS.search(combined) and METAL.search(combined))
    quant = bool(QUANT.search(combined))
    if REVIEW_TITLE.search(title):
        category, reason = 'preliminary_exclusion', 'review_or_meta_analysis_explicit_in_title'
    elif not title or not rec['abstract'].strip():
        category, reason = 'manual_review', 'missing_title_or_abstract'
    elif (bio or sorp) and quant:
        category, reason = 'candidate', 'phenotype_and_quantitative_terms_in_title_or_abstract'
    else:
        category, reason = 'manual_review', 'insufficient_evidence_for_automated_classification'
    return dict(rec, screen_category=category, screen_reason=reason,
                biofilm_signal=int(bio), biosorption_signal=int(sorp),
                quantitative_signal=int(quant), screen_method='title_abstract_rules_v1')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='data/raw/pubmed_records.tsv')
    parser.add_argument('--outdir', default='data/interim')
    parser.add_argument('--reportdir', default='results/reports')
    args = parser.parse_args()
    infile = Path(args.input)
    if not infile.is_file():
        parser.error(f'Input not found: {infile}')
    outdir, reportdir = Path(args.outdir), Path(args.reportdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reportdir.mkdir(parents=True, exist_ok=True)
    seen = set()
    rows = []
    with infile.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or any(col not in reader.fieldnames for col in FIELDS):
            parser.error('Input columns do not match 01_pubmed_search.py output')
        for n, rec in enumerate(reader, start=2):
            pmid = rec['pmid'].strip()
            if not pmid or not pmid.isdigit() or pmid in seen:
                parser.error(f'Missing, invalid or duplicate PMID on input line {n}: {pmid!r}')
            seen.add(pmid)
            rows.append(screen({field: rec[field] or '' for field in FIELDS}))
    output = outdir / 'pubmed_screening.tsv'
    with output.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS + EXTRA, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(r['screen_category'] for r in rows)
    summary = reportdir / 'pubmed_screening_summary.tsv'
    with summary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['category','records'])
        for category in ['candidate','manual_review','preliminary_exclusion']:
            writer.writerow([category, counts[category]])
        writer.writerow(['TOTAL', len(rows)])
    print(f'Screened {len(rows)} records. Original file unchanged.')
    print(f'Screening: {output}')
    print(f'Summary: {summary}')
    print('IMPORTANT: preliminary exclusions and candidates require human validation.')

if __name__ == '__main__':
    main()
