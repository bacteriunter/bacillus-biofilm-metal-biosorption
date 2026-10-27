#!/usr/bin/env python3
"""Conservative, read-only triage of existing phenotype audit and locally saved full texts.
No automatic phenotype verification; no network, no genome reanalysis.
"""
import argparse, csv, json, re, time, urllib.request, urllib.parse, urllib.error, xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

MEASUREMENT = re.compile(r'\b(?:biofilm (?:formation|production|biomass|assay|quantification)|crystal violet|microtiter|microtitre|optical density|OD\s*570|OD\s*595|exopolysaccharide (?:production|yield)|EPS (?:production|yield)|biosorption capacity|adsorption capacity|sorption capacity|metal (?:removal|immobilization)|Langmuir|q\s*max)\b', re.I)
VALUE = re.compile(r'(?<!\w)\d+(?:\.\d+)?\s*(?:mg\s*/\s*g|mg\s*g\s*[-−]1|%|OD\s*\d+|µg\s*/\s*mL|μg\s*/\s*mL)', re.I)
NEGATIVE = re.compile(r'\b(?:inhibit(?:ed|ion|s)?|anti[- ]biofilm|against|pathogen biofilm|biofilm of (?:Pseudomonas|Staphylococcus|Candida)|suppressed biofilm)\b', re.I)
CONDITIONS = re.compile(r'\b(?:NaCl|salt tolerance|thermotolerance|hemolytic|antibiotic susceptibility|sequence identity)\b', re.I)

def norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower())

def snippets(text, strain, max_items=4):
    # Require a strain mention within a 450-character window around a measurement keyword.
    found=[]
    aliases={strain, strain.replace('-', ''), strain.replace('-', ' ')}
    for m in MEASUREMENT.finditer(text):
        start=max(0,m.start()-220); end=min(len(text),m.end()+230)
        seg=' '.join(text[start:end].split())
        if not any(re.search(r'(?<![A-Za-z0-9])'+re.escape(a)+r'(?![A-Za-z0-9])',seg,re.I) for a in aliases if len(a)>=3):
            continue
        found.append((seg, bool(VALUE.search(seg)), bool(NEGATIVE.search(seg))))
        if len(found)>=max_items: break
    return found


def retrieve_oa(pmid, destination):
    """Retrieve only Europe PMC open-access full-text XML, no publisher scraping."""
    base='https://www.ebi.ac.uk/europepmc/webservices/rest/'
    def get(url):
        req=urllib.request.Request(url, headers={'User-Agent':'BMetalPhenotypeAudit/1.0 (academic research)', 'Accept':'application/json, application/xml'})
        with urllib.request.urlopen(req, timeout=25) as response:
            return response.read()
    try:
        q=urllib.parse.urlencode({'query':f'EXT_ID:{pmid} AND SRC:MED','format':'json','pageSize':1})
        meta=json.loads(get(base+'search?'+q))
        hits=meta.get('resultList',{}).get('result',[])
        if not hits: return 'metadata_not_found',''
        hit=hits[0]
        pmcid=hit.get('pmcid','')
        if not pmcid or str(hit.get('isOpenAccess','')).upper()!='Y':
            return 'not_open_access',pmcid
        raw=get(base+urllib.parse.quote(pmcid)+'/fullTextXML')
        doc=ET.fromstring(raw)
        # Preserve paragraphs, titles and individual table cells; exclude references.
        for back in list(doc.iter('back')):
            back.clear()
        fragments=[]
        for el in doc.iter():
            if el.tag in ('title','p','td','th'):
                content=' '.join(' '.join(el.itertext()).split())
                if content: fragments.append(content)
        if not fragments: return 'fulltext_xml_without_text',pmcid
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text('\n'.join(fragments)+'\n',encoding='utf-8')
        return 'retrieved_oa',pmcid
    except urllib.error.HTTPError as e:
        return f'http_{e.code}',''
    except (urllib.error.URLError,TimeoutError,ET.ParseError,ValueError) as e:
        return 'retrieval_error_'+type(e).__name__,''

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--top',type=int,default=20)
    ap.add_argument('--fetch-oa-top',type=int,default=0,help='Fetch OA full texts for top N unique PMIDs from Europe PMC')
    args=ap.parse_args()
    root=Path(args.root).resolve()
    source=root/'results/functional/phenotype_audit/phenotype_audit.tsv'
    if not source.is_file():
        ap.error(f'No existe {source}')
    with source.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    # Publication/strain is the review unit. Do not collapse distinct assembly candidates silently.
    grouped=defaultdict(list)
    for r in rows:
        pmid=(r.get('pmid') or '').strip()
        strain=(r.get('strain') or '').strip()
        if pmid and strain: grouped[(pmid,r.get('species','').strip(),strain)].append(r)
    out=[]
    fulltext_dir=root/'data/raw/biosorption_fulltext'
    oa_dir=root/'data/raw/phenotype_fulltext'
    for (pmid,species,strain),rs in grouped.items():
        domains=sorted({d.strip() for r in rs for d in (r.get('domains') or '').split(';') if d.strip()})
        assemblies=sorted({r.get('assembly_accession','').strip() for r in rs if r.get('assembly_accession','').strip()})
        excerpts=list(dict.fromkeys(r.get('evidence_excerpt','') for r in rs if r.get('evidence_excerpt','')))
        local=[]
        for folder in (fulltext_dir,oa_dir):
            for ext in ('.txt','.xml'):
                path=folder/(pmid+ext)
                if path.is_file():
                    local.append(path)
        evidence=[]
        for path in local:
            content=path.read_text(encoding='utf-8',errors='replace')
            for seg,has_value,negative in snippets(content,strain):
                evidence.append((str(path.relative_to(root)),seg,has_value,negative))
            if evidence: break # avoid duplicate txt+xml
        excerpt=' | '.join(excerpts)[:1800]
        text=evidence[0][1] if evidence else excerpt
        has_endpoint=bool(MEASUREMENT.search(text))
        has_value=bool(VALUE.search(text))
        negative=bool(NEGATIVE.search(text))
        conditions=bool(CONDITIONS.search(text))
        # A ranking for human review, never a verification or phenotype assignment.
        score=(4 if evidence else 0)+(2 if has_endpoint else 0)+(2 if has_value else 0)+(1 if len(assemblies)==1 else 0)-(3 if negative else 0)-(1 if conditions else 0)
        status='local_fulltext_context_found' if evidence else ('local_fulltext_no_strain_endpoint_context' if local else 'no_local_fulltext')
        out.append(dict(pmid=pmid,species=species,strain=strain,assemblies=';'.join(assemblies),domains=';'.join(domains),
                        candidate_assembly_count=len(assemblies),local_fulltext=';'.join(str(p.relative_to(root)) for p in local),
                        text_status=status,priority_score=score,measurement_term_in_excerpt_or_context=int(has_endpoint),
                        numeric_value_in_excerpt_or_context=int(has_value),possible_other_organism_or_inhibition=int(negative),
                        context_excerpt=text[:800],manual_decision='pending_primary_source_review',
                        phenotype_value_verified='',phenotype_unit_verified='',phenotype_method_verified='',
                        strain_genome_identity_verified='',reviewer_notes=''))
    out.sort(key=lambda r:(-r['priority_score'],r['pmid'],r['strain']))
    retrieval=[]
    if args.fetch_oa_top>0:
        chosen=[]; seen=set()
        for r in out:
            if r['pmid'] not in seen:
                seen.add(r['pmid']); chosen.append(r['pmid'])
            if len(chosen)>=args.fetch_oa_top: break
        for pmid in chosen:
            target=oa_dir/(pmid+'.txt')
            if target.is_file() or (fulltext_dir/(pmid+'.txt')).is_file() or (fulltext_dir/(pmid+'.xml')).is_file():
                status,pmcid='already_local',''
            else:
                status,pmcid=retrieve_oa(pmid,target)
                time.sleep(0.4)
            retrieval.append({'pmid':pmid,'pmcid':pmcid,'status':status})
            print(f'PMID {pmid}: {status} {pmcid}',flush=True)
        # Re-score all units whose publication was downloaded.
        for r in out:
            pmid=r['pmid']; strain=r['strain']
            local=[p for p in (oa_dir/(pmid+'.txt'),fulltext_dir/(pmid+'.txt'),fulltext_dir/(pmid+'.xml')) if p.is_file()]
            if not local: continue
            r['local_fulltext']=';'.join(str(p.relative_to(root)) for p in local)
            evidence=[]
            for path in local:
                evidence=snippets(path.read_text(encoding='utf-8',errors='replace'),strain)
                if evidence: break
            if evidence:
                seg,has_value,negative=evidence[0]
                r['context_excerpt']=seg[:800]
                r['text_status']='local_fulltext_context_found'
                r['measurement_term_in_excerpt_or_context']=1
                r['numeric_value_in_excerpt_or_context']=int(has_value)
                r['possible_other_organism_or_inhibition']=int(negative)
                r['priority_score']=4+2+2*int(has_value)+int(r['candidate_assembly_count']==1)-3*int(negative)-int(bool(CONDITIONS.search(seg)))
            else:
                r['text_status']='local_fulltext_no_strain_endpoint_context'
        out.sort(key=lambda r:(-r['priority_score'],r['pmid'],r['strain']))

    dest=root/'results/functional/phenotype_audit'
    dest.mkdir(parents=True,exist_ok=True)
    path=dest/'targeted_review_shortlist.tsv'
    cols=list(out[0]) if out else []
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(out)
    if retrieval:
        with (dest/'oa_retrieval_manifest.tsv').open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['pmid','pmcid','status'],delimiter='\t'); w.writeheader(); w.writerows(retrieval)
    report=dict(oa_retrieval={'attempted_unique_pmids':len(retrieval),'statuses':{status:sum(x['status']==status for x in retrieval) for status in sorted({x['status'] for x in retrieval})}},generated_utc=datetime.now(timezone.utc).isoformat(),input_audit_rows=len(rows),unique_publication_strain_units=len(out),
                with_local_fulltext=sum(bool(x['local_fulltext']) for x in out),
                with_strain_and_endpoint_context_in_local_fulltext=sum(x['text_status']=='local_fulltext_context_found' for x in out),
                verified_new_phenotypes=0,notice='Ranking is heuristic, not phenotype validation; manually inspect primary source tables, figures and methods.')
    (dest/'targeted_review_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))
    print('\nPRIMEROS CANDIDATOS (revisión manual, no verificados):')
    for r in out[:args.top]:
        print(f"PMID {r['pmid']} | {r['species']} {r['strain']} | puntuación {r['priority_score']} | {r['text_status']} | {r['assemblies']}")
    print('\nArchivos:',path, dest/'targeted_review_summary.json')

if __name__=='__main__': main()
