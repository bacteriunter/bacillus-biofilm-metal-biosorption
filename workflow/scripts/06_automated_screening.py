#!/usr/bin/env python3
"""Exploratory, reproducible PubMed title/abstract screening; no automatic eligibility decisions."""
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'data/raw/pubmed_records.tsv'
OUT = ROOT / 'data/interim/pubmed_automated_evidence.tsv'
REPORT = ROOT / 'results/reports/pubmed_automated_evidence_summary.json'

PATTERNS = {
    'bacillus': r'\bbacillus\b|\bb\.\s*(?:subtilis|licheniformis|velezensis|cereus|megaterium|safensis|altitudinis|amyloliquefaciens)\b',
    'biofilm': r'\bbiofilms?\b|\bbiofilm[- ]forming\b',
    'matrix_eps': r'\bextracellular (?:polymeric substances?|polysaccharides?|matrix|polymers?)\b|\bexopolysaccharides?\b|\bEPS\b|\bextracellular matrix\b',
    'biosorption': r'\bbiosorp\w*\b|\bbioadsorp\w*\b|\b(?:metal|ion|cadmium|lead|copper|nickel|chromium|zinc|mercury) (?:uptake|adsorption|binding|sorption|removal)\b',
    'metal': r'\b(?:heavy metals?|cadmium|lead|chromium|copper|nickel|zinc|mercury|arsenic|manganese|strontium|Pb\(II\)|Cd\(II\)|Cr\(VI\)|Cu\(II\)|Ni\(II\))\b',
    'plant': r'\b(?:phytoremediation|rhizosph\w*|plant[- ]assisted|plant growth|seedlings?|roots?)\b',
    'quantitative': r'\b(?:quantif\w*|measur\w*|assay|kinetic\w*|isotherm\w*|capacity|efficien\w*|mg/g|mg L[-−]?1|percent|percentage|removal rate|crystal violet)\b',
    'resistance': r'\b(?:metal resistance|metal tolerance|resistant to|tolerant to)\b',
    'isolated_eps': r'\b(?:purified|isolated|extracted)\s+(?:extracellular|EPS|exopolysaccharide|polymeric)\b',
}
COMPILED = {k: re.compile(v, re.I) for k, v in PATTERNS.items()}
SPECIES = re.compile(r'\bBacillus\s+(?:sp\.?|[a-z][a-z-]+)(?:\s+(?:strain\s+)?[A-Z0-9][A-Za-z0-9_.-]{1,18})?', re.I)
REVIEW = re.compile(r'\b(?:systematic review|scoping review|meta-analysis|a review|review of|review on)\b', re.I)


def main():
    if not SOURCE.is_file():
        raise SystemExit(f'Missing input: {SOURCE}')
    with SOURCE.open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or not {'pmid', 'title', 'abstract'}.issubset(rows[0]):
        raise SystemExit('Input missing records or required columns')
    seen = set()
    counts = Counter()
    output = []
    for r in rows:
        pmid = (r.get('pmid') or '').strip()
        if not pmid or pmid in seen:
            raise SystemExit(f'Duplicate or blank PMID: {pmid!r}')
        seen.add(pmid)
        title = r.get('title') or ''
        abstract = r.get('abstract') or ''
        combined = title + '\n' + abstract
        hits = {k: bool(p.search(combined)) for k, p in COMPILED.items()}
        possible_review = bool(REVIEW.search(title))
        if not abstract.strip():
            tier = 'insufficient_metadata'
        elif not hits['bacillus']:
            tier = 'context_or_uncertain'
        elif (hits['biofilm'] or hits['matrix_eps']) and hits['biosorption'] and hits['metal']:
            tier = 'cross_domain_candidate'
        elif hits['biofilm'] or hits['matrix_eps'] or (hits['biosorption'] and hits['metal']):
            tier = 'single_domain_candidate'
        else:
            tier = 'context_or_uncertain'
        if possible_review:
            tier = 'possible_review_check'
        counts[tier] += 1
        mentions = list(dict.fromkeys(m.group(0) for m in SPECIES.finditer(combined)))
        output.append({**r,
            'evidence_tier': tier,
            'matched_signals': ';'.join(k for k, v in hits.items() if v),
            'taxon_mentions_unverified': ';'.join(mentions[:15]),
            'possible_review_title': int(possible_review),
            'eligibility_status': 'not_assessed',
            'method_version': 'title_abstract_evidence_v1',
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(output[0]), delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(output)
    report = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source': str(SOURCE.relative_to(ROOT)),
        'records': len(rows),
        'evidence_tiers': dict(sorted(counts.items())),
        'automatic_exclusions': 0,
        'eligibility_decisions': 0,
        'method': 'fixed transparent regex evidence flags on title and abstract',
        'limitations': 'Not a validated eligibility classifier; taxon mentions are unverified and no strain-genome linkage is inferred.',
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'Processed {len(rows)} records. No eligibility exclusions.')
    print(f'Evidence table: {OUT.relative_to(ROOT)}')
    print(f'Summary: {REPORT.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
