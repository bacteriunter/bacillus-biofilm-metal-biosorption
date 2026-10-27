#!/usr/bin/env python3
"""Supplement missing gene families from UniProtKB without modifying module 26 outputs.
Exact gene symbol matches are checked against UniProtKB gene-name fields.
"""
import argparse, csv, hashlib, json, re, time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

MISSING = {
    'epsA':'eps_matrix','epsB':'eps_matrix',
    'pgsB':'poly_gamma_glutamate','pgsC':'poly_gamma_glutamate',
    'pgsE':'poly_gamma_glutamate','capE':'poly_gamma_glutamate',
    'zntA':'metal_homeostasis_not_biosorption','zitB':'metal_homeostasis_not_biosorption'
}
BASE='https://rest.uniprot.org/uniprotkb/search'
FIELDS=['family','mechanism','gene','protein_accession','reference_organism','evidence_note','source_url','length_aa','sequence_sha256','review_status','query','sequence_version']

def request_json(url):
    for attempt in range(4):
        try:
            req=Request(url,headers={'User-Agent':'bmetal-reference-panel/1.0 (academic research)'})
            with urlopen(req,timeout=60) as response:
                return json.load(response)
        except Exception:
            if attempt == 3: raise
            time.sleep(2**attempt)

def exact_genes(entry):
    found=set()
    for g in entry.get('genes',[]):
        for k in ('geneName','orderedLocusNames','orfNames','synonyms'):
            x=g.get(k,[])
            if isinstance(x,dict): x=[x]
            for v in x:
                if isinstance(v,dict) and v.get('value'): found.add(v['value'].lower())
    return found

def fetch(query):
    # Query only B. subtilis 168 first; second search uses reviewed Bacillus entries.
    url=BASE+'?'+urlencode({'query':query,'format':'json','size':500})
    seen=set()
    while url:
        if url in seen: raise RuntimeError('Repeated pagination URL')
        seen.add(url)
        payload=request_json(url)
        for entry in payload.get('results',[]): yield entry
        # UniProt JSON search provides next-page URL in Link HTTP header, not in JSON.
        # Requests are deliberately scoped to a single page; report truncation if 500.
        if len(payload.get('results',[]))==500:
            raise RuntimeError(f'Query may exceed 500 records; refine query: {query}')
        url=None

def seq_valid(s):
    return bool(s) and re.fullmatch('[ACDEFGHIKLMNPQRSTVWYBXZJUO]+',s) is not None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--panel',default='data/reference/functional_panel')
    args=ap.parse_args()
    root=Path(args.panel); root.mkdir(parents=True,exist_ok=True)
    existing=root/'reference_catalog.tsv'
    if not existing.exists(): raise SystemExit('Missing original reference_catalog.tsv')
    with existing.open(newline='') as f: original=list(csv.DictReader(f,delimiter='\t'))
    if not original: raise SystemExit('Original catalog empty')
    supplement=[]; log=[]; seen_acc={r.get('protein_accession','') for r in original}
    for gene,mechanism in MISSING.items():
        # Unreviewed sequences only for exact gene names in strain 168; avoid unsupervised genus-wide TrEMBL.
        queries=[f'gene_exact:{gene} AND organism_id:224308',
                 f'gene_exact:{gene} AND taxonomy_id:1386 AND reviewed:true']
        count=0
        for q in queries:
            entries=list(fetch(q))
            accepted=0
            for e in entries:
                acc=e.get('primaryAccession',''); seq=e.get('sequence',{}).get('value','').upper()
                if gene.lower() not in exact_genes(e) or not seq_valid(seq) or not acc or acc in seen_acc: continue
                org=e.get('organism',{}).get('scientificName','')
                if not org.startswith('Bacillus '): continue
                seen_acc.add(acc)
                row={'family':mechanism,'mechanism':mechanism,'gene':gene,'protein_accession':acc,
                     'reference_organism':org,'evidence_note':'Exact gene symbol; function not experimentally validated by this retrieval',
                     'source_url':f'https://rest.uniprot.org/uniprotkb/{acc}',
                     'length_aa':len(seq),'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                     'review_status':e.get('entryType',''),'query':q,
                     'sequence_version':e.get('entryAudit',{}).get('sequenceVersion','')}
                supplement.append((row,seq)); accepted+=1
            log.append({'gene':gene,'query':q,'returned':len(entries),'accepted':accepted})
            count+=accepted
        print(f'{gene}\t{count}',flush=True)
    out=root/'reference_supplement_26b.tsv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',extrasaction='ignore');w.writeheader()
        for r,_ in supplement:w.writerow(r)
    with (root/'reference_supplement_26b.faa').open('w') as f:
        for r,seq in supplement:
            f.write(f">{r['protein_accession']}|{r['gene']}\n")
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
    (root/'reference_supplement_26b_provenance.json').write_text(json.dumps({'queries':log,'note':'Gene symbol match is not orthology or experimentally proven function; existing panel untouched'},indent=2)+'\n')
    print(f'ADDED: {len(supplement)} sequences; original files unchanged')
    print('NEXT: inspect supplemental catalog and integrate with original panel in the workflow')

if __name__=='__main__':main()
