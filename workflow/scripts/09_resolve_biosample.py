#!/usr/bin/env python3
"""Retrieve BioSample attributes for provisional Assembly hits; no automatic verification."""
import argparse, csv, json, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime, timezone

INPUT = Path('data/interim/ncbi_assembly_candidates.tsv')
OUTPUT = Path('data/interim/ncbi_biosample_resolution.tsv')
CACHE = Path('data/interim/ncbi_biosample_metadata_cache.jsonl')
SUMMARY = Path('results/reports/ncbi_biosample_resolution_summary.json')

def request_xml(accession, email):
    args = urllib.parse.urlencode({'db':'biosample','id':accession,'retmode':'xml','tool':'bmetal','email':email})
    req = urllib.request.Request('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+args, headers={'User-Agent':'bmetal/0.1 (BioSample metadata research)'})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=35) as response:
                return response.read()
        except Exception:
            if attempt == 3: raise
            time.sleep(2 ** attempt)

def parse_biosample(payload):
    root = ET.fromstring(payload)
    node = root if root.tag == 'BioSample' else root.find('.//BioSample')
    if node is None: raise ValueError('BioSample node absent')
    attrs = {}
    for attr in node.findall('.//Attribute'):
        key = (attr.attrib.get('harmonized_name') or attr.attrib.get('attribute_name') or '').strip().lower()
        if key: attrs.setdefault(key, []).append((attr.text or '').strip())
    return {'accession':node.attrib.get('accession',''), 'attributes':attrs,
            'organism': (node.find('.//Organism').attrib.get('taxonomy_name','') if node.find('.//Organism') is not None else '')}

def normalized(value):
    return re.sub(r'\s+', ' ', value.strip().casefold())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--email', required=True)
    args = parser.parse_args()
    if not INPUT.exists(): parser.error(f'Missing {INPUT}')
    with INPUT.open(newline='', encoding='utf-8') as fh:
        rows = [r for r in csv.DictReader(fh, delimiter='\t') if r['match_status'] == 'ambiguous' and r['biosample_accession'].strip()]
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    existing = {}
    if CACHE.exists():
        for line in CACHE.read_text(encoding='utf-8').splitlines():
            try:
                obj = json.loads(line)
                if obj.get('biosample_accession') and obj.get('status') == 'ok': existing[obj['biosample_accession']] = obj
            except json.JSONDecodeError: pass
    accessions = sorted(set(r['biosample_accession'].strip() for r in rows))
    with CACHE.open('a', encoding='utf-8') as cachefile:
        for i, accession in enumerate(accessions, 1):
            if accession in existing: continue
            try:
                parsed = parse_biosample(request_xml(accession, args.email))
                obj = {'biosample_accession':accession,'status':'ok', **parsed}
            except Exception as exc:
                obj = {'biosample_accession':accession,'status':'error','error':str(exc)}
            cachefile.write(json.dumps(obj, ensure_ascii=False)+'\n'); cachefile.flush()
            existing[accession] = obj
            if i % 25 == 0: print(f'Processed {i}/{len(accessions)} BioSamples', flush=True)
            time.sleep(0.37)
    fields = ['species_mention','strain_candidate','pmids','dois','assembly_accession','biosample_accession','biosample_organism','biosample_strain','biosample_isolate','biosample_culture_collection','resolution_status','resolution_basis']
    counts = {}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t'); writer.writeheader()
        for row in rows:
            obj = existing[row['biosample_accession'].strip()]
            attrs = obj.get('attributes', {})
            strain = attrs.get('strain', []); isolate = attrs.get('isolate', []); collection = attrs.get('culture_collection', [])
            values = [v for v in strain + isolate + collection if v]
            candidate = normalized(row['strain_candidate'])
            matching = [v for v in values if normalized(v) == candidate]
            if obj['status'] != 'ok': status, basis = 'unresolved', 'biosample_fetch_error'
            elif matching: status, basis = 'metadata_match_candidate', 'exact_attribute_value'
            elif values: status, basis = 'unresolved', 'different_or_nonmatching_identifiers_not_proof_of_conflict'
            else: status, basis = 'unresolved', 'no_strain_identifiers_in_biosample'
            counts[status] = counts.get(status, 0) + 1
            writer.writerow({'species_mention':row['species_mention'],'strain_candidate':row['strain_candidate'],'pmids':row['pmids'],'dois':row['dois'],'assembly_accession':row['assembly_accession'],'biosample_accession':row['biosample_accession'],'biosample_organism':obj.get('organism',''),'biosample_strain':'; '.join(strain),'biosample_isolate':'; '.join(isolate),'biosample_culture_collection':'; '.join(collection),'resolution_status':status,'resolution_basis':basis})
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),'input':str(INPUT),'rows_processed':len(rows),'unique_biosamples':len(accessions),'resolution_status_rows':counts,'output':str(OUTPUT),'cache':str(CACHE),'limitations':'Exact attribute matches are candidates only; strain codes can be reused. Different identifiers are not automatically conflicts. No phenotype-genome link verified.'}, indent=2)+'\n',encoding='utf-8')
    print(f'Wrote {OUTPUT} and {SUMMARY}')

if __name__ == '__main__': main()
