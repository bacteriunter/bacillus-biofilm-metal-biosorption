#!/usr/bin/env python3
"""Prepare a reproducible, overlap-stratified manual screening sample.
No eligibility decisions are inferred or made.
"""
import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = ['pmid','doi','title','abstract','search_groups','manual_decision','exclusion_reason','reviewer_notes']

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',default='data/raw/pubmed_records.tsv')
    p.add_argument('--output',default='data/interim/manual_validation_sample.tsv')
    p.add_argument('--report',default='results/reports/manual_validation_sample.json')
    p.add_argument('--n',type=int,default=300)
    p.add_argument('--seed',type=int,default=20261023)
    a=p.parse_args()
    if a.n<1: p.error('--n must be positive')
    src=Path(a.input)
    with src.open(encoding='utf-8',newline='') as f:
        records=list(csv.DictReader(f,delimiter='\t'))
    required={'pmid','doi','title','abstract','search_groups'}
    if not records or not required.issubset(records[0]):
        raise SystemExit('Missing records or required columns')
    ids=[r['pmid'] for r in records]
    if len(ids)!=len(set(ids)):
        raise SystemExit('Duplicate PMIDs detected; stop and audit input')
    if a.n>len(records): p.error('--n exceeds population')
    strata=defaultdict(list)
    for r in records:
        groups=tuple(sorted(set(g.strip() for g in r['search_groups'].replace('|',';').split(';') if g.strip())))
        strata['+'.join(groups) or 'unassigned'].append(r)
    # Proportional allocation using largest remainder, exact total n.
    raw={k:a.n*len(v)/len(records) for k,v in strata.items()}
    alloc={k:int(raw[k]) for k in strata}
    remaining=a.n-sum(alloc.values())
    for k in sorted(strata,key=lambda k:(-(raw[k]-alloc[k]),k))[:remaining]: alloc[k]+=1
    rng=random.Random(a.seed)
    sampled=[]
    for k in sorted(strata):
        sampled.extend(rng.sample(sorted(strata[k],key=lambda r:r['pmid']),alloc[k]))
    rng.shuffle(sampled)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists(): raise SystemExit(f'Output exists, refusing overwrite: {out}')
    with out.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader()
        for r in sampled:
            w.writerow({**{k:r.get(k,'') for k in FIELDS},'manual_decision':'','exclusion_reason':'','reviewer_notes':''})
    report=Path(a.report);report.parent.mkdir(parents=True,exist_ok=True)
    if report.exists(): raise SystemExit(f'Report exists, refusing overwrite: {report}')
    data={'population':len(records),'sample_size':len(sampled),'seed':a.seed,
          'strata_population':dict(sorted((k,len(v)) for k,v in strata.items())),
          'strata_sample':dict(sorted(alloc.items())),
          'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
          'sampling_method':'proportional stratified random sampling by exact search-group combination; largest remainder allocation',
          'eligibility_decisions':0,'note':'Sample for human screening; not an exclusion set or validated classifier.'}
    report.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f'Created {out} ({len(sampled)} records)')
    print(f'Report: {report}')
    print('No automatic eligibility decisions made.')
if __name__=='__main__': main()
