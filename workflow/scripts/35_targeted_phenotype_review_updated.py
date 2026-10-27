#!/usr/bin/env python3
"""Conservative, read-only triage of existing phenotype audit and locally saved full texts.
No automatic phenotype verification; no network, no genome reanalysis.
"""
import argparse, csv, json, re
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

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--top',type=int,default=20)
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
    for (pmid,species,strain),rs in grouped.items():
        domains=sorted({d.strip() for r in rs for d in (r.get('domains') or '').split(';') if d.strip()})
        assemblies=sorted({r.get('assembly_accession','').strip() for r in rs if r.get('assembly_accession','').strip()})
        excerpts=list(dict.fromkeys(r.get('evidence_excerpt','') for r in rs if r.get('evidence_excerpt','')))
        local=[]
        for ext in ('.txt','.xml'):
            path=fulltext_dir/(pmid+ext)
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
    dest=root/'results/functional/phenotype_audit'
    dest.mkdir(parents=True,exist_ok=True)
    path=dest/'targeted_review_shortlist.tsv'
    cols=list(out[0]) if out else []
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t'); w.writeheader(); w.writerows(out)
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),input_audit_rows=len(rows),unique_publication_strain_units=len(out),
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
