#!/usr/bin/env python3
"""Non-exclusionary, auditable prioritization of PubMed screening records."""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

PATTERNS = {
    'biofilm': r'\bbiofilms?\b|\bpellicles?\b|\bsessile\b',
    'matrix': r'extracellular matrix|extracellular polymeric|\bEPS\b|\bexopolysaccharides?\b|\btasA\b|\bepsA',
    'biosorption': r'biosorp\w*|bioadsorp\w*|metal (?:binding|uptake|adsorption)|sorption (?:of |by )?(?:heavy )?metals?',
    'metal': r'\b(?:cadmium|lead|chromium|copper|nickel|zinc|arsenic|mercury|cobalt|uranium|\bCd\b|\bPb\b|\bCr\b|\bCu\b|\bNi\b|\bZn\b|\bAs\b|\bHg\b)\b',
    'quantitative': r'\b(?:quantif\w*|measur\w*|assay\w*|capacity|isotherm\w*|kinetic\w*|removal efficien\w*|crystal violet|adsorption capacit\w*|\bmg/g\b|\bmg l-1\b)\b',
    'plant': r'\b(?:plant\w*|rhizosph\w*|phytoremediation|root\w*|soil\w*)\b',
    'strain': r'\bstrain\b|\bisolate\b|\b[A-Z]{1,4}[- ]?\d{2,}\b',
}
WEIGHTS = {'biofilm':2,'matrix':3,'biosorption':3,'metal':2,'quantitative':2,'plant':1,'strain':1}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='data/interim/pubmed_screening.tsv')
    ap.add_argument('--output',default='data/interim/pubmed_prioritized.tsv')
    ap.add_argument('--summary',default='results/reports/pubmed_prioritization_summary.tsv')
    args=ap.parse_args()
    src=Path(args.input)
    if not src.exists(): ap.error(f'Missing input: {src}')
    out=Path(args.output); summary=Path(args.summary)
    out.parent.mkdir(parents=True,exist_ok=True);summary.parent.mkdir(parents=True,exist_ok=True)
    counts=Counter(); total=0
    with src.open(encoding='utf-8-sig',newline='') as fin, out.open('w',encoding='utf-8',newline='') as fout:
        reader=csv.DictReader(fin,delimiter='\t')
        if not reader.fieldnames: ap.error('Input has no header')
        required={'pmid','title','abstract','screen_category'}
        missing=required-set(reader.fieldnames)
        if missing: ap.error('Missing columns: '+', '.join(sorted(missing)))
        extras=['priority_score','priority_tier','priority_signals','priority_method']
        writer=csv.DictWriter(fout,fieldnames=reader.fieldnames+extras,delimiter='\t',extrasaction='ignore')
        writer.writeheader()
        for row in reader:
            text=' '.join([row.get('title') or '',row.get('abstract') or ''])
            flags={name:bool(re.search(pattern,text,re.IGNORECASE)) for name,pattern in PATTERNS.items()}
            score=sum(WEIGHTS[name] for name,flag in flags.items() if flag)
            # Priority is a review order, never an eligibility decision.
            if flags['biosorption'] and flags['metal'] and (flags['matrix'] or flags['biofilm']): tier='high'
            elif (flags['biofilm'] and flags['matrix']) or (flags['biosorption'] and flags['metal']): tier='medium'
            else: tier='low'
            row.update(priority_score=str(score),priority_tier=tier,
                       priority_signals=';'.join(k for k,v in flags.items() if v),
                       priority_method='keyword_prioritization_v1_nonexclusionary')
            writer.writerow(row);counts[tier]+=1;total+=1
    with summary.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['priority_tier','records'])
        for tier in ('high','medium','low'):w.writerow([tier,counts[tier]])
        w.writerow(['TOTAL',total])
    print(f'Prioritized {total} records; no records excluded.')
    print(f'Output: {out}\nSummary: {summary}')

if __name__=='__main__': main()
