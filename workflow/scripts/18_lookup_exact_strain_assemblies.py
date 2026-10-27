#!/usr/bin/env python3
"""Look up NCBI Assembly candidates for exact strain identifiers via E-utilities.
No candidate is marked confirmed without explicit matching strain metadata.
"""
import csv, json, time, urllib.parse, urllib.request
from pathlib import Path

OUT=Path('data/processed/strain_assembly_lookup.tsv')
REPORT=Path('results/reports/strain_assembly_lookup_summary.json')
QUERIES=[('Bacillus salmalaya','139SI','ATCC BAA-2268'),('Bacillus subtilis','EG-Pb5',''),('Bacillus','MRS-1','ATCC 55673')]
FIELDS=['requested_species','requested_strain','alias','assembly_accession','biosample_accession','ncbi_organism','ncbi_infraspecies','assembly_title','match_status','ncbi_uid']
BASE='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(endpoint,params):
    params.update({'tool':'bmetal_strain_lookup','email':'research@example.org'})
    url=BASE+endpoint+'?'+urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url,timeout=30) as r:return json.load(r)
        except Exception:
            if attempt==2:raise
            time.sleep(2*(attempt+1))

def find_candidates(strain,alias):
    terms=[f'"{strain}"[All Fields]']
    if alias:terms.append(f'"{alias}"[All Fields]')
    uids=set()
    for term in terms:
        data=get('esearch.fcgi',{'db':'assembly','term':term,'retmode':'json','retmax':'100'})
        uids.update(data.get('esearchresult',{}).get('idlist',[]))
        time.sleep(.4)
    return sorted(uids)

def normalize(s):return ''.join(c.lower() for c in str(s) if c.isalnum())

def main():
    OUT.parent.mkdir(parents=True,exist_ok=True);REPORT.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    for species,strain,alias in QUERIES:
        print(f'Consultando {species} {strain} ...',flush=True)
        try:
            uids=find_candidates(strain,alias)
            if not uids:
                rows.append(dict(zip(FIELDS,[species,strain,alias,'','','','','','not_found',''])))
                continue
            for i in range(0,len(uids),20):
                data=get('esummary.fcgi',{'db':'assembly','id':','.join(uids[i:i+20]),'retmode':'json','report':'full'})
                result=data.get('result',{})
                for uid in uids[i:i+20]:
                    r=result.get(uid,{})
                    accession=r.get('assemblyaccession','')
                    biosample=r.get('biosampleaccn','')
                    organism=r.get('organism','')
                    infra=r.get('infraspecieslist',[]) or []
                    infra_text='; '.join(str(x) for x in infra)
                    title=r.get('assemblyname','')
                    blob=' '.join([str(organism),infra_text,str(title)])
                    matches=normalize(strain) in normalize(blob) or (alias and normalize(alias) in normalize(blob))
                    status='exact_identifier_in_metadata' if matches else 'search_hit_unverified'
                    rows.append(dict(zip(FIELDS,[species,strain,alias,accession,biosample,organism,infra_text,title,status,uid])))
                time.sleep(.4)
        except Exception as e:
            rows.append(dict(zip(FIELDS,[species,strain,alias,'','','','','','lookup_error: '+str(e),''])))
    with OUT.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)
    counts={}
    for row in rows:counts[row['match_status']]=counts.get(row['match_status'],0)+1
    REPORT.write_text(json.dumps({'queries':len(QUERIES),'candidate_rows':len(rows),'status_counts':counts,'output':str(OUT),'note':'Identifier in Assembly metadata is only a candidate; verify BioSample identity before assigning phenotype.'},indent=2),encoding='utf-8')
    print(json.dumps(counts,indent=2))
    print('Resultados:',OUT)
if __name__=='__main__':main()
