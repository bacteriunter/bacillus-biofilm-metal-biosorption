#!/usr/bin/env python3
"""Extract auditable phenotype evidence from PubMed titles/abstracts.

Text mining only: numerical mentions are NOT validated phenotype measurements.
No automatic eligibility decisions or strain-to-measurement associations.
"""
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'data/raw/pubmed_records.tsv'
OUT = ROOT / 'data/interim/pubmed_phenotype_evidence.tsv'
NUMOUT = ROOT / 'data/interim/pubmed_quantitative_mentions.tsv'
REPORT = ROOT / 'results/reports/pubmed_phenotype_extraction_summary.json'

DOMAINS = {
    'biofilm': r'\bbiofilms?\b|\bbiofilm[- ]form(?:ation|ing)\b|\bcrystal violet\b|\badherent biomass\b',
    'eps_matrix': r'\bEPS\b|\bextracellular polymeric substances?\b|\bextracellular polysaccharides?\b|\bexopolysaccharides?\b|\bextracellular matrix\b',
    'metal_sorption': r'\bbiosorp\w*\b|\bbioadsorp\w*\b|\badsorp\w*\b|\bsorp\w*\b|\bmetal[- ]bind\w*\b|\bmetal uptake\b|\bmetal removal\b',
    'plant_context': r'\bphytoremediat\w*\b|\brhizosph\w*\b|\bplant[- ]assisted\b|\bplant growth\b|\broot\b',
}
METALS = {
    'Cd': r'\bcadmium\b|\bCd\s*\(?II\)?\b|\bCd2\+',
    'Pb': r'\blead\b|\bPb\s*\(?II\)?\b|\bPb2\+',
    'Cr': r'\bchromium\b|\bchromate\b|\bCr\s*\(?VI\)?\b|\bCr\s*\(?III\)?\b',
    'Cu': r'\bcopper\b|\bCu\s*\(?II\)?\b',
    'Ni': r'\bnickel\b|\bNi\s*\(?II\)?\b',
    'Zn': r'\bzinc\b|\bZn\s*\(?II\)?\b',
    'Hg': r'\bmercury\b|\bHg\s*\(?II\)?\b',
    'As': r'\barsenic\b|\barsenate\b|\barsenite\b',
    'Mn': r'\bmanganese\b|\bMn\s*\(?II\)?\b',
}
# Detect candidate numerical mentions, retaining surrounding text; never interpret as measured endpoints automatically.
QUANT = re.compile(
    r'(?<![\w.])(?P<value>\d+(?:[.,]\d+)?)\s*'
    r'(?P<unit>%|mg\s*/\s*g|mg\s*g\s*[-−]?1|mg\s*/\s*L|mg\s*L\s*[-−]?1|'
    r'g\s*/\s*L|g\s*L\s*[-−]?1|µg\s*/\s*mL|μg\s*/\s*mL|'
    r'OD\s*\d{3}|nm|h|hours?|days?|min|mM|µM|μM)(?!\w)', re.IGNORECASE)
SPLIT = re.compile(r'(?<=[.!?])\s+(?=[A-Z(])')

def hits(text, patterns):
    return [key for key, pattern in patterns.items() if re.search(pattern, text, re.I)]

def main():
    if not SOURCE.exists():
        raise SystemExit(f'Missing input: {SOURCE}')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    evidence_fields = ['pmid','doi','year','search_groups','source_field','sentence_index','evidence_text','domains','metals','has_numeric_mention','interpretation_status']
    num_fields = ['pmid','doi','source_field','sentence_index','value_text','unit_text','context_text','domains','metals','measurement_status']
    stats = Counter()
    seen_pmids = set()
    with SOURCE.open(encoding='utf-8-sig', newline='') as fi, OUT.open('w', encoding='utf-8', newline='') as fo, NUMOUT.open('w', encoding='utf-8', newline='') as fn:
        reader = csv.DictReader(fi, delimiter='\t')
        required = {'pmid','doi','year','title','abstract','search_groups'}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f'Missing columns: {sorted(required - set(reader.fieldnames or []))}')
        ew = csv.DictWriter(fo, fieldnames=evidence_fields, delimiter='\t'); ew.writeheader()
        nw = csv.DictWriter(fn, fieldnames=num_fields, delimiter='\t'); nw.writeheader()
        for row in reader:
            stats['records'] += 1
            pmid = row['pmid'].strip()
            if pmid in seen_pmids:
                stats['duplicate_pmid_rows'] += 1
            seen_pmids.add(pmid)
            record_domains = set()
            for source_field in ('title','abstract'):
                for idx, sentence in enumerate(SPLIT.split((row.get(source_field) or '').strip()), 1):
                    if not sentence: continue
                    domains = hits(sentence, DOMAINS)
                    if not domains: continue
                    record_domains.update(domains)
                    metals = hits(sentence, METALS)
                    quantities = list(QUANT.finditer(sentence))
                    ew.writerow({'pmid':pmid,'doi':row['doi'],'year':row['year'], 'search_groups':row['search_groups'],
                                 'source_field':source_field,'sentence_index':idx,'evidence_text':sentence,
                                 'domains':';'.join(domains),'metals':';'.join(metals),
                                 'has_numeric_mention':int(bool(quantities)), 'interpretation_status':'unverified_text_evidence'})
                    stats['evidence_rows'] += 1
                    for match in quantities:
                        nw.writerow({'pmid':pmid,'doi':row['doi'],'source_field':source_field,'sentence_index':idx,
                                     'value_text':match.group('value'),'unit_text':match.group('unit'),
                                     'context_text':sentence,'domains':';'.join(domains),'metals':';'.join(metals),
                                     'measurement_status':'unverified_numeric_mention'})
                        stats['numeric_mention_rows'] += 1
            if record_domains:
                stats['records_with_evidence'] += 1
                for domain in record_domains: stats[f'records_{domain}'] += 1
    report = {'generated_utc':datetime.now(timezone.utc).isoformat(), 'source':str(SOURCE.relative_to(ROOT)),
              **dict(sorted(stats.items())), 'evidence_output':str(OUT.relative_to(ROOT)),
              'numeric_output':str(NUMOUT.relative_to(ROOT)),
              'limitations':'Title/abstract regex extraction only. Numeric mentions may be conditions, times, doses or unrelated results. No strain-to-value linkage, eligibility decisions or comparable phenotypes established.'}
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(f"Processed {stats['records']} records; {stats['evidence_rows']} evidence passages; {stats['numeric_mention_rows']} numeric mentions.")
    print(f'Report: {REPORT.relative_to(ROOT)}')

if __name__ == '__main__': main()
