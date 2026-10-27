#!/usr/bin/env python3
"""Reproducible UniProt reviewed-reference retrieval; no genome annotation.
Run from project root. Python standard library only.
"""
import argparse,csv,hashlib,json,re,time,urllib.parse,urllib.request
from pathlib import Path

# Explicit candidate scope; reference evidence is recorded, not presumed.
PANEL={
 'eps_matrix':['epsA','epsB','epsC','epsD','epsE','epsF','epsG','epsH','epsI','epsJ','epsK','epsL','epsM','epsN','epsO'],
 'protein_matrix':['tapA','sipW','tasA','bslA'],
 'biofilm_regulation':['sinR','sinI','spo0A','abrB','degU','comA','slrR','remA','remB'],
 'poly_gamma_glutamate':['pgsA','pgsB','pgsC','pgsE','capA','capB','capC','capE'],
 'cell_wall_polymers':['tagO','tagA','tagB','tagD','tagF','dltA','dltB','dltC','dltD'],
 'metal_homeostasis_not_biosorption':['czcD','cadA','copA','copZ','csoR','arsB','arsC','arsR','mntH','mntR','zntA','znuA','znuB','znuC','zitB','fur'],
}
FIELDS='accession,id,reviewed,protein_name,gene_primary,organism_name,organism_id,sequence,protein_existence'
API='https://rest.uniprot.org/uniprotkb/search'

def fetch(query,timeout,retries):
    params=urllib.parse.urlencode({'query':query,'format':'tsv','fields':FIELDS,'size':'500'})
    url=API+'?'+params
    for attempt in range(retries):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'bmetal-research-reference-panel/1.0','Accept':'text/tab-separated-values'})
            with urllib.request.urlopen(req,timeout=timeout) as resp:
                body=resp.read().decode('utf-8')
                link=resp.headers.get('Link','')
            yield body,url
            while link:
                m=re.search(r'<([^>]+)>;\s*rel="next"',link)
                if not m: break
                url=m.group(1)
                req=urllib.request.Request(url,headers={'User-Agent':'bmetal-research-reference-panel/1.0'})
                with urllib.request.urlopen(req,timeout=timeout) as resp:
                    body=resp.read().decode('utf-8');link=resp.headers.get('Link','')
                yield body,url
            return
        except Exception:
            if attempt+1==retries: raise
            time.sleep(2**attempt)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',default='data/reference/functional_panel')
    ap.add_argument('--timeout',type=int,default=90)
    ap.add_argument('--retries',type=int,default=3)
    args=ap.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    records={}; searches=[]
    for group,genes in PANEL.items():
        for gene in genes:
            # Reviewed proteins within genus Bacillus. No arbitrary per-gene count limit.
            query=f'(taxonomy_id:1386) AND (reviewed:true) AND (gene_exact:{gene})'
            found=0
            for body,url in fetch(query,args.timeout,args.retries):
                reader=csv.DictReader(body.splitlines(),delimiter='\t')
                for row in reader:
                    acc=row.get('Entry','').strip();seq=row.get('Sequence','').strip().upper()
                    if not acc or not seq: continue
                    # Restrict to exact primary gene where present; keep other hits separately out of scope.
                    primary=row.get('Gene Names (primary)','').strip()
                    if primary and gene.lower() not in [g.lower() for g in primary.split()]: continue
                    key=(group,gene,acc)
                    records[key]={
                      'group':group,'query_gene':gene,'uniprot_accession':acc,
                      'entry_name':row.get('Entry Name',''),
                      'protein_name':row.get('Protein names',''),
                      'primary_gene':primary,
                      'organism':row.get('Organism',''),
                      'taxon_id':row.get('Organism (ID)',''),
                      'protein_existence':row.get('Protein existence',''),
                      'sequence_length':len(seq),
                      'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
                      'source_url':'https://rest.uniprot.org/uniprotkb/'+acc,
                      'sequence':seq}
                    found+=1
            searches.append({'group':group,'gene':gene,'query':query,'reviewed_records':found})
            print(f'{group}\t{gene}\t{found}',flush=True)
            time.sleep(0.15)
    if not records: raise SystemExit('ERROR: no reviewed reference sequences found; outputs not overwritten')
    cols=[c for c in next(iter(records.values())) if c!='sequence']
    with (out/'reference_catalog.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
        for k in sorted(records): w.writerow({c:records[k][c] for c in cols})
    with (out/'reference_proteins.faa').open('w') as f:
        for k in sorted(records):
            r=records[k]
            f.write(f">{r['uniprot_accession']}|{r['group']}|{r['query_gene']}\n")
            seq=r['sequence']
            for i in range(0,len(seq),80):f.write(seq[i:i+80]+'\n')
    (out/'reference_retrieval.json').write_text(json.dumps({'source':'UniProtKB reviewed','taxon':'Bacillus taxon 1386','method':'gene_exact with primary-gene verification','searches':searches,'records':len(records),'limitations':['Reviewed-only coverage is incomplete for many Bacillus lineages','Gene name matching is not functional verification or orthology','Metal homeostasis proteins are not biosorption markers','This is a traceable seed reference library, not a biologically complete curated panel']},indent=2)+'\n')
    print(f'READY: {len(records)} proteins; {out}/reference_proteins.faa')
if __name__=='__main__':main()
