#!/usr/bin/env python3
"""Extract traceable biosorption-capacity candidates from available full texts.
No candidate is automatically verified or attributed to a strain.
Standard library only. Run from the bmetal project root.
"""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

CAPACITY = re.compile(r'(?<![\w.])(?P<value>\d+(?:[.,]\d+)?)\s*(?P<unit>mg\s*(?:[A-Za-z0-9()+\-]{0,15}\s*)?[/⁄]\s*g(?:\s*[-−]?1)?|mg\s*g\s*[-−]1|mg\s+per\s+g)(?!\w)', re.I)
BIOSORPTION = re.compile(r'biosorp|adsorp|sorp|uptake|binding|bound|capacity|q\s*max|qmax|langmuir|remov', re.I)
METAL = re.compile(r'\b(?:Cd|Pb|Cu|Zn|Ni|Cr|As|Hg|Co|Mn|Fe|Au|Ag)(?:\s*\(?\s*(?:II|III|IV|VI|[2-6]\+?)\s*\)?)?\b|cadmium|lead|copper|zinc|nickel|chromium|arsenic|mercury|cobalt|manganese|iron|gold|silver', re.I)
DRY = re.compile(r'dry\s*(?:weight|mass|biomass|cell)|dried|lyophiliz|freeze.dri', re.I)
FRESH = re.compile(r'fresh\s*(?:weight|mass|biomass|cell)|wet\s*(?:weight|mass|biomass)', re.I)
EPS = re.compile(r'\bEPS\b|extracellular polymer|exopolysaccharide|poly.gamma.glutam|\bPGA\b', re.I)
SPORE = re.compile(r'\bspores?\b', re.I)
DOSE = re.compile(r'initial|dosage|dose|ratio|concentration|added|input|loading', re.I)
MAXIMUM = re.compile(r'maximum|\bmax\b|q\s*max|langmuir|capacity', re.I)
SENTENCE = re.compile(r'(?<=[.!?])\s+(?=[A-Z0-9(])')
FIELDS = ['pmid','doi','source_file','source_type','sentence_index','value','unit_original','context','metal_mention','capacity_language','dose_language','biomass_basis','material_flag','species_strain_candidates','assembly_candidates','biosample_candidates','evidence_status','review_decision','verified_strain','verified_assembly','verified_capacity_mg_g_dry','verification_notes']

def read_tsv(path):
    if not path.exists(): return []
    with path.open(encoding='utf-8-sig',newline='') as fh:
        return list(csv.DictReader(fh,delimiter='\t'))

def read_text(path):
    raw=path.read_text(encoding='utf-8',errors='replace')
    if path.suffix.lower() in ('.xml','.nxml'):
        try:
            root=ET.fromstring(raw)
            for node in root.iter():
                if node.tag.rsplit('}',1)[-1] in ('ref-list','table-wrap','fig','supplementary-material'):
                    node.clear()
            raw=' '.join(root.itertext())
        except ET.ParseError: return ''
    elif path.suffix.lower()=='.json':
        try:
            obj=json.loads(raw)
            if isinstance(obj,dict):
                raw=' '.join(str(obj.get(k,'')) for k in ('title','abstract','body','fulltext','text','content'))
        except json.JSONDecodeError: return ''
    return re.sub(r'\s+',' ',raw).strip()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.',help='Project root')
    args=ap.parse_args()
    root=Path(args.root)
    fulltext=root/'data/raw/biosorption_fulltext'
    if not fulltext.is_dir(): raise SystemExit(f'Missing full-text directory: {fulltext}')
    lit=read_tsv(root/'data/raw/pubmed_records.tsv')
    links=read_tsv(root/'data/processed/genome_phenotype_publication_links.tsv')
    mentions=read_tsv(root/'data/interim/bacillus_taxon_mentions.tsv')
    bypmid=defaultdict(list)
    bydoi={}
    for r in lit:
        p=str(r.get('pmid','')).strip(); d=str(r.get('doi','')).strip().lower()
        if p: bypmid[p].append(r)
        if d and p: bydoi[d]=p
    ml=defaultdict(set)
    for r in mentions:
        p=str(r.get('pmid','')).strip()
        if p and r.get('strain_candidate'):
            ml[p].add((r.get('species_mention','')+' '+r.get('strain_candidate','')).strip())
    al=defaultdict(set); bl=defaultdict(set)
    for r in links:
        p=str(r.get('pmid','')).strip()
        if not p: continue
        for key, dest in [('assembly_accession',al),('biosample_accession',bl)]:
            v=str(r.get(key,'')).strip()
            if v: dest[p].add(v)
    manifest=read_tsv(root/'data/interim/biosorption_fulltext_manifest.tsv')
    file_to_pmid={}
    for row in manifest:
        pmid=str(row.get('pmid','')).strip()
        if not pmid: continue
        for key in ('pmcid','file','filepath','file_path','output','output_path','local_path','path','filename'):
            value=str(row.get(key,'')).strip()
            if value:
                file_to_pmid[value.lower()]=pmid
                file_to_pmid[Path(value).stem.lower()]=pmid
                file_to_pmid[Path(value).name.lower()]=pmid
    files=sorted(p for p in fulltext.rglob('*') if p.is_file() and p.suffix.lower() in ('.xml','.nxml','.txt','.json'))
    out=[]; unmapped=[]; stats=Counter()
    for file in files:
        text=read_text(file)
        if not text: stats['unreadable_or_empty']+=1; continue
        pmids=set(re.findall(r'(?<!\d)(\d{7,9})(?!\d)',file.stem)) & set(bypmid)
        if not pmids:
            # Only use exact PMID from the file name, never assume a PMID from arbitrary article text.
            pmids={p for p in bypmid if p in file.stem and not re.search(r'\d'+re.escape(p)+r'\d',file.stem)}
        if not pmids:
            for key in (file.stem.lower(),file.name.lower(),str(file.relative_to(root)).lower()):
                if key in file_to_pmid: pmids.add(file_to_pmid[key])
        if not pmids:
            # PMC filenames may use PMCID; retain for separate manifest reconciliation.
            unmapped.append(str(file.relative_to(root))); stats['files_without_pmid_mapping']+=1
            continue
        for pmid in sorted(pmids):
            doi=next((r.get('doi','') for r in bypmid[pmid]),'')
            sentences=SENTENCE.split(text)
            for i,sentence in enumerate(sentences):
                if not BIOSORPTION.search(sentence): continue
                for match in CAPACITY.finditer(sentence):
                    context=' '.join(sentences[max(0,i-1):min(len(sentences),i+2)])[:1500]
                    # Never assert that a mg/g number is a sorption capacity: context is for review.
                    basis='fresh_or_wet' if FRESH.search(context) else 'dry_or_lyophilized' if DRY.search(context) else 'unspecified'
                    material='eps_or_polymer' if EPS.search(context) else 'spore' if SPORE.search(context) else 'not_specified'
                    out.append(dict(pmid=pmid,doi=doi,source_file=str(file.relative_to(root)),source_type='fulltext',sentence_index=i,value=match.group('value').replace(',','.'),unit_original=match.group('unit'),context=context,metal_mention='; '.join(sorted(set(METAL.findall(context)),key=str.lower)),capacity_language='yes' if MAXIMUM.search(context) else 'no',dose_language='yes' if DOSE.search(context) else 'no',biomass_basis=basis,material_flag=material,species_strain_candidates='; '.join(sorted(ml[pmid])),assembly_candidates='; '.join(sorted(al[pmid])),biosample_candidates='; '.join(sorted(bl[pmid])),evidence_status='unverified_fulltext_candidate',review_decision='pending',verified_strain='',verified_assembly='',verified_capacity_mg_g_dry='',verification_notes=''))
            stats['files_mapped']+=1
    dest=root/'data/processed'; dest.mkdir(parents=True,exist_ok=True)
    output=dest/'fulltext_biosorption_capacity_candidates.tsv'
    with output.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t'); w.writeheader(); w.writerows(out)
    report=root/'results/reports/fulltext_biosorption_capacity_summary.json'; report.parent.mkdir(parents=True,exist_ok=True)
    summary={'generated_utc':datetime.now(timezone.utc).isoformat(),'fulltext_files_seen':len(files),'fulltext_files_mapped':stats['files_mapped'],'unmapped_files':unmapped,'candidate_mentions':len(out),'unique_pmids_with_candidates':len({r['pmid'] for r in out}),'verified_pairs':0,'output':str(output.relative_to(root)),'limitations':'Full text heuristic extraction; no capacity, strain or genome pairing automatically validated. Unmapped PMCID files must be reconciled using the retrieval manifest; no articles excluded.'}
    report.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
