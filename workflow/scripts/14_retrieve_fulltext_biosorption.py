#!/usr/bin/env python3
"""Retrieve openly available full texts for the 53 candidate biosorption papers.

Uses Europe PMC REST; preserves provenance. This is a retrieval step, not
phenotype validation or strain-genome linkage. Requires internet in WSL.
"""
import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INFILE = ROOT / 'data/processed/biosorption_evidence_candidates.tsv'
OUTDIR = ROOT / 'data/raw/biosorption_fulltext'
MANIFEST = ROOT / 'data/interim/biosorption_fulltext_manifest.tsv'
REPORT = ROOT / 'results/reports/biosorption_fulltext_retrieval_summary.json'
BASE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/'
HEADERS = {'User-Agent': 'bmetal-research/1.0 (Europe PMC public REST; academic use)'}


def fetch(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=35) as resp:
        return resp.read()


def get_metadata(pmid):
    url = BASE + 'search?query=' + urllib.parse.quote(f'EXT_ID:{pmid} AND SRC:MED') + '&format=json&pageSize=2'
    payload = json.loads(fetch(url))
    results = payload.get('resultList', {}).get('result', [])
    return next((r for r in results if str(r.get('id')) == str(pmid)), {})


def get_text(xmlbytes):
    root = ET.fromstring(xmlbytes)
    # Extract article paragraphs, captions, and tables with section provenance.
    records = []
    for el in root.iter():
        if el.tag.rsplit('}', 1)[-1] not in {'p', 'table-wrap', 'caption'}:
            continue
        content = ' '.join(''.join(el.itertext()).split())
        if len(content) >= 30:
            records.append(content)
    return '\n'.join(dict.fromkeys(records))


def main():
    if not INFILE.exists():
        raise SystemExit(f'Missing: {INFILE}')
    with INFILE.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not {'pmid', 'doi'}.issubset(reader.fieldnames or []):
            raise SystemExit('Input must contain pmid and doi')
        studies = {}
        for row in reader:
            pmid = row['pmid'].strip()
            if pmid:
                studies.setdefault(pmid, row['doi'].strip())
    OUTDIR.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (pmid, doi) in enumerate(sorted(studies.items()), 1):
        xmlpath = OUTDIR / f'{pmid}.xml'
        txtpath = OUTDIR / f'{pmid}.txt'
        status = 'unavailable'
        pmcid = ''
        error = ''
        try:
            metadata = get_metadata(pmid)
            pmcid = metadata.get('pmcid', '') or ''
            if not metadata:
                status = 'metadata_not_found'
            elif not pmcid:
                status = 'no_pmcid'
            elif xmlpath.exists() and txtpath.exists():
                status = 'cached_fulltext'
            else:
                url = BASE + urllib.parse.quote(pmcid) + '/fullTextXML'
                try:
                    xml = fetch(url)
                    text = get_text(xml)
                    if not text:
                        status = 'empty_fulltext'
                    else:
                        xmlpath.write_bytes(xml)
                        txtpath.write_text(text + '\n', encoding='utf-8')
                        status = 'retrieved_fulltext'
                except urllib.error.HTTPError as exc:
                    status = 'fulltext_not_available' if exc.code in (401, 403, 404) else 'http_error'
                    error = f'HTTP {exc.code}'
        except (urllib.error.URLError, TimeoutError, ValueError, ET.ParseError) as exc:
            status = 'retrieval_error'
            error = str(exc)[:180]
        rows.append({'pmid': pmid, 'doi': doi, 'pmcid': pmcid,
                     'status': status, 'fulltext_xml': str(xmlpath.relative_to(ROOT)) if xmlpath.exists() else '',
                     'fulltext_txt': str(txtpath.relative_to(ROOT)) if txtpath.exists() else '',
                     'error': error})
        print(f'[{i}/{len(studies)}] PMID {pmid}: {status}', flush=True)
        time.sleep(0.35)
    fields = ['pmid','doi','pmcid','status','fulltext_xml','fulltext_txt','error']
    with MANIFEST.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    report = {'generated_utc': datetime.now(timezone.utc).isoformat(),
              'publications_attempted': len(studies),
              'statuses': dict(Counter(row['status'] for row in rows)),
              'manifest': str(MANIFEST.relative_to(ROOT)),
              'note': 'Openly accessible full texts only; retrieval does not establish strain-phenotype-genome links.'}
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('Report:', REPORT.relative_to(ROOT))

if __name__ == '__main__':
    main()
