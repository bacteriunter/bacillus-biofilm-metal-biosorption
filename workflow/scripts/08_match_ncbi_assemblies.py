#!/usr/bin/env python3
"""Conservative NCBI Assembly candidate matching for literature-derived Bacillus strains.

Uses E-utilities; results are discovery candidates, NOT validated genotype–phenotype links.
Standard library only. Cache enables safe resumption. No genome downloads.
"""
import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
import urllib.error
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
FIELDS = ['pmid','doi','year','title','search_groups','mention_text','species_mention','strain_candidate','mention_level','source_field','verification_status']
OUTFIELDS = ['species_mention','strain_candidate','pmids','dois','query','assembly_accession','assembly_uid','ncbi_species','ncbi_infraspecies','biosample_accession','match_status','match_basis']

def norm(s):
    return re.sub(r'[^a-z0-9]+', '', (s or '').lower())

def get_json(endpoint, params, email, delay, retries=4):
    params.update({'tool':'bmetal_literature_genome_linkage','email':email,'retmode':'json'})
    url = BASE + endpoint + '?' + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            time.sleep(delay)
            req = urllib.request.Request(url, headers={'User-Agent':f'bmetal/0.1 ({email})'})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if attempt == retries - 1: raise RuntimeError(f'NCBI request failed: {e}') from e
            time.sleep(2 ** attempt)

def ncbi_search(species, strain, email, delay):
    # Search retrieves possible matches; a field match in a search result is NOT proof of strain identity.
    term = f'"{species}"[Organism] AND "{strain}"[All Fields]'
    result = get_json('esearch.fcgi', {'db':'assembly','term':term,'retmax':'20'}, email, delay)
    ids = result.get('esearchresult', {}).get('idlist', [])
    count = int(result.get('esearchresult', {}).get('count', 0))
    if not ids: return term, count, []
    summary = get_json('esummary.fcgi', {'db':'assembly','id':','.join(ids),'report':'full'}, email, delay)
    block = summary.get('result', {})
    return term, count, [block.get(str(uid), {}) for uid in ids]

def infra_text(entry):
    v = entry.get('infraspecieslist', [])
    if isinstance(v, list):
        return '; '.join(f"{x.get('sub_type','')}={x.get('sub_value','')}" for x in v if isinstance(x, dict))
    return str(v or '')

def evaluate(entry, species, strain):
    ncbi_species = str(entry.get('speciesname') or entry.get('organism') or '')
    infra = infra_text(entry)
    # Only an explicit strain subtype is a strong candidate; still requires BioSample verification.
    strain_values = []
    for item in entry.get('infraspecieslist', []) if isinstance(entry.get('infraspecieslist', []), list) else []:
        if isinstance(item, dict) and str(item.get('sub_type','')).lower() == 'strain':
            strain_values.append(str(item.get('sub_value','')))
    species_ok = norm(species) == norm(ncbi_species) or norm(species) in norm(ncbi_species)
    exact = species_ok and any(norm(strain) == norm(v) for v in strain_values)
    return ('exact_candidate' if exact else 'ambiguous',
            'assembly_infraspecies_strain_equal' if exact else 'search_hit_not_verified', infra, ncbi_species)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--email', required=True, help='Real contact email required by NCBI')
    p.add_argument('--input', default='data/interim/bacillus_taxon_mentions.tsv')
    p.add_argument('--output', default='data/interim/ncbi_assembly_candidates.tsv')
    p.add_argument('--report', default='results/reports/ncbi_assembly_candidates_summary.json')
    p.add_argument('--cache', default='data/interim/ncbi_assembly_query_cache.jsonl')
    p.add_argument('--delay', type=float, default=0.38, help='Minimum delay between NCBI requests without API key')
    args = p.parse_args()
    if '@' not in args.email: p.error('Provide a valid contact email with --email')
    if args.delay < 0.34: p.error('Use --delay >= 0.34 without an NCBI API key')
    with open(args.input, encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = set(FIELDS) - set(reader.fieldnames or [])
        if missing: p.error('Missing columns: ' + ', '.join(sorted(missing)))
        rows = list(reader)
    pairs = {}
    skipped = Counter()
    for r in rows:
        if r['mention_level'].strip() != 'candidate_strain':
            skipped['not_candidate_strain'] += 1; continue
        species = re.sub(r'\s+', ' ', r['species_mention'].strip())
        strain = re.sub(r'\s+', ' ', r['strain_candidate'].strip())
        if not re.fullmatch(r'Bacillus\s+[A-Za-z][A-Za-z-]+', species, re.I):
            skipped['unresolved_species'] += 1; continue
        if not (2 <= len(strain) <= 35 and re.search(r'[A-Za-z]', strain) and re.search(r'\d', strain)):
            skipped['weak_strain_code'] += 1; continue
        key = (species.lower(), strain.lower())
        pairs.setdefault(key, {'species':species,'strain':strain,'pmids':set(),'dois':set()})
        if r['pmid']: pairs[key]['pmids'].add(r['pmid'])
        if r['doi']: pairs[key]['dois'].add(r['doi'])
    cache_path = Path(args.cache); cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache = {}
    if cache_path.exists():
        with cache_path.open(encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    obj = json.loads(line)
                    cache[obj['key']] = obj
    output_rows = []
    for index, (key, pair) in enumerate(sorted(pairs.items()), 1):
        cache_key = '|'.join(key)
        if cache_key in cache:
            record = cache[cache_key]
        else:
            term, count, entries = ncbi_search(pair['species'], pair['strain'], args.email, args.delay)
            record = {'key':cache_key,'query':term,'total_hits':count,'entries':entries,
                      'retrieved_utc':datetime.now(timezone.utc).isoformat()}
            with cache_path.open('a', encoding='utf-8') as f:
                f.write(json.dumps(record, ensure_ascii=False) + '\n')
            cache[cache_key] = record
        entries = record['entries']
        if not entries:
            entries = [None]
        for entry in entries:
            status, basis, infra, ncbi_species = evaluate(entry, pair['species'], pair['strain']) if entry else ('not_found','no_assembly_search_hits','','')
            output_rows.append({
                'species_mention':pair['species'], 'strain_candidate':pair['strain'],
                'pmids':';'.join(sorted(pair['pmids'])), 'dois':';'.join(sorted(pair['dois'])),
                'query':record['query'], 'assembly_accession':(entry or {}).get('assemblyaccession',''),
                'assembly_uid':(entry or {}).get('uid',''), 'ncbi_species':ncbi_species,
                'ncbi_infraspecies':infra, 'biosample_accession':(entry or {}).get('biosampleaccn',''),
                'match_status':status, 'match_basis':basis})
        if index % 25 == 0 or index == len(pairs):
            print(f'Processed {index}/{len(pairs)} unique species–strain queries', flush=True)
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=OUTFIELDS, delimiter='\t'); w.writeheader(); w.writerows(output_rows)
    report = {'source':args.input,'literature_rows':len(rows),'unique_species_strain_queries':len(pairs),
              'candidate_assembly_rows':len(output_rows),'match_status_rows':dict(Counter(r['match_status'] for r in output_rows)),
              'skipped_mentions':dict(skipped),'cache':args.cache,
              'note':'Search results are provisional. Even exact_candidate is NOT a verified strain–phenotype linkage. Validate BioSample, paper methods and genome metadata before modeling.'}
    report_path = Path(args.report); report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(f'Candidates: {output}\nSummary: {report_path}')

if __name__ == '__main__': main()
