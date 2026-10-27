#!/usr/bin/env python3
"""Reproducible PubMed retrieval for bmetal (NCBI E-utilities)."""
import argparse
import csv
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
QUERIES = {
    'biofilm': '(Bacillus[Title/Abstract]) AND (biofilm*[Title/Abstract] OR "extracellular matrix"[Title/Abstract] OR "extracellular polymeric substances"[Title/Abstract])',
    'biosorption': '(Bacillus[Title/Abstract]) AND (biosorp*[Title/Abstract] OR bioadsorp*[Title/Abstract] OR "metal adsorption"[Title/Abstract] OR "metal binding"[Title/Abstract] OR "metal uptake"[Title/Abstract])',
    'plant_metal': '(Bacillus[Title/Abstract]) AND (phytoremediat*[Title/Abstract] OR "plant-assisted remediation"[Title/Abstract] OR "plant growth-promoting"[Title/Abstract]) AND (metal*[Title/Abstract] OR cadmium[Title/Abstract] OR lead[Title/Abstract] OR chromium[Title/Abstract] OR copper[Title/Abstract] OR nickel[Title/Abstract] OR zinc[Title/Abstract])',
}
FIELDS = ['pmid', 'doi', 'title', 'authors', 'year', 'journal', 'abstract', 'search_groups']

def request(endpoint, params, email, tool, api_key=None):
    params = {**params, 'tool': tool, 'email': email}
    if api_key:
        params['api_key'] = api_key
    url = BASE + endpoint + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': f'{tool}/1.0 ({email})'}), timeout=90) as response:
                payload = response.read()
            time.sleep(0.36 if not api_key else 0.12)
            return payload
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 4:
                raise RuntimeError(f'NCBI request failed: {endpoint}') from exc
            time.sleep(min(2 ** attempt, 16))

def plain(node):
    return ' '.join(''.join(node.itertext()).split()) if node is not None else ''

def article_record(article):
    citation = article.find('./MedlineCitation')
    pmid = plain(citation.find('./PMID'))
    pub = citation.find('./Article')
    ids = article.findall('./PubmedData/ArticleIdList/ArticleId')
    doi = next((plain(item) for item in ids if item.get('IdType') == 'doi'), '')
    title = plain(pub.find('./ArticleTitle'))
    abstract = ' '.join(plain(x) for x in pub.findall('./Abstract/AbstractText'))
    authors = []
    for a in pub.findall('./AuthorList/Author'):
        collective = plain(a.find('./CollectiveName'))
        name = ' '.join(x for x in (plain(a.find('./ForeName')), plain(a.find('./LastName'))) if x)
        if collective or name:
            authors.append(collective or name)
    year = plain(pub.find('./Journal/JournalIssue/PubDate/Year'))
    if not year:
        medline_date = plain(pub.find('./Journal/JournalIssue/PubDate/MedlineDate'))
        year = medline_date[:4] if medline_date[:4].isdigit() else ''
    if not year:
        year = plain(pub.find('./ArticleDate/Year'))
    return {'pmid': pmid, 'doi': doi, 'title': title, 'authors': '; '.join(authors),
            'year': year, 'journal': plain(pub.find('./Journal/Title')),
            'abstract': abstract, 'search_groups': ''}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--email', required=True, help='Your contact email required by NCBI')
    parser.add_argument('--outdir', default='data/raw')
    parser.add_argument('--reportdir', default='results/reports')
    parser.add_argument('--api-key', default=os.getenv('NCBI_API_KEY'))
    parser.add_argument('--batch-size', type=int, default=200)
    args = parser.parse_args()
    if args.batch_size < 1 or args.batch_size > 500:
        parser.error('--batch-size must be between 1 and 500')
    outdir, reportdir = Path(args.outdir), Path(args.reportdir)
    outdir.mkdir(parents=True, exist_ok=True)
    reportdir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    records, group_ids, log = {}, {}, {'retrieved_at_utc': started, 'database': 'PubMed', 'queries': {}, 'ncbi_api': 'E-utilities'}
    for group, query in QUERIES.items():
        root = ET.fromstring(request('esearch.fcgi', {'db': 'pubmed', 'term': query, 'retmode': 'xml', 'retmax': 0, 'usehistory': 'y'}, args.email, 'bmetal_pubmed', args.api_key))
        count = int(plain(root.find('./Count')))
        webenv, query_key = plain(root.find('./WebEnv')), plain(root.find('./QueryKey'))
        if count and (not webenv or not query_key):
            raise RuntimeError(f'Missing PubMed history token for {group}')
        group_ids[group] = set()
        log['queries'][group] = {'query': query, 'pubmed_count': count}
        print(f'{group}: {count} PubMed records', flush=True)
        for start in range(0, count, args.batch_size):
            root = ET.fromstring(request('efetch.fcgi', {'db': 'pubmed', 'query_key': query_key, 'WebEnv': webenv,
                'retstart': start, 'retmax': min(args.batch_size, count-start), 'retmode': 'xml'}, args.email, 'bmetal_pubmed', args.api_key))
            for article in root.findall('./PubmedArticle'):
                rec = article_record(article)
                if rec['pmid']:
                    group_ids[group].add(rec['pmid'])
                    records.setdefault(rec['pmid'], rec)
            print(f'  retrieved batch starting {start}', flush=True)
        log['queries'][group]['unique_pmids_retrieved'] = len(group_ids[group])
    for pmid, rec in records.items():
        rec['search_groups'] = ';'.join(group for group in QUERIES if pmid in group_ids[group])
    log['unique_pmids_total'] = len(records)
    log['completed_at_utc'] = datetime.now(timezone.utc).isoformat()
    with (outdir / 'pubmed_records.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(records[pmid] for pmid in sorted(records, key=lambda p: int(p)))
    with (outdir / 'pubmed_search_log.json').open('w', encoding='utf-8') as f:
        json.dump(log, f, ensure_ascii=False, indent=2)
    with (reportdir / 'pubmed_search_summary.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['search_group', 'pubmed_count', 'unique_pmids_retrieved'])
        for group in QUERIES:
            entry = log['queries'][group]
            writer.writerow([group, entry['pubmed_count'], entry['unique_pmids_retrieved']])
        writer.writerow(['TOTAL_UNIQUE', '', len(records)])
    print(f'Done: {len(records)} unique records. Output: {outdir / "pubmed_records.tsv"}')

if __name__ == '__main__':
    main()
