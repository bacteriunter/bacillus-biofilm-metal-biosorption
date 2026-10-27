#!/usr/bin/env python3
"""Extract candidate Bacillus taxon/strain mentions from PubMed records.

Exploratory text mining only: mentions are not verified experimental strains.
Uses Python standard library. Never overwrites the input file.
"""
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SOURCE = Path('data/raw/pubmed_records.tsv')
OUT = Path('data/interim/bacillus_taxon_mentions.tsv')
SUMMARY = Path('results/reports/bacillus_taxon_mentions_summary.json')

# Conservative strain-designator extraction. A designation is only a candidate.
SPECIES = re.compile(r'\bBacillus\s+(?P<species>[a-z][a-z-]{2,}|sp\.)\b', re.I)
ABBREV = re.compile(r'\bB\.\s*(?P<species>[a-z][a-z-]{2,})\b')
STRAIN = re.compile(
    r'\s+(?:(?:strain|isolate|subsp\.?|str\.?)\s+)?'
    r'(?P<strain>(?=[A-Za-z0-9_-]*\d)[A-Za-z][A-Za-z0-9_-]{1,18}|'
    r'[A-Z]{2,8}[-_]?[0-9]{1,8})\b'
)
STOP = {'and', 'with', 'from', 'under', 'using', 'isolated', 'biofilm', 'metal',
        'biosorption', 'producing', 'cells', 'strain', 'strains', 'sp', 'spp'}
FIELDS = ['pmid', 'doi', 'year', 'title', 'search_groups', 'mention_text',
          'species_mention', 'strain_candidate', 'mention_level', 'source_field',
          'verification_status']


def mentions(text, source_field):
    hits = []
    for pattern, abbreviated in ((SPECIES, False), (ABBREV, True)):
        for match in pattern.finditer(text):
            species = match.group('species').lower()
            if species in STOP or len(species) < 3:
                continue
            suffix = text[match.end():match.end() + 35]
            # Explicit 'strain/isolate' or an immediately following alphanumeric code.
            explicit = re.match(r'\s+(?:strain|isolate|str\.)\s+([A-Za-z0-9][\w-]{1,18})\b', suffix, re.I)
            code = explicit.group(1) if explicit else None
            if not code:
                m = STRAIN.match(suffix)
                if m:
                    code = m.group('strain')
            if code and code.lower() in STOP:
                code = None
            # Do not assign designations to generic 'Bacillus sp.' unless explicitly indicated.
            if species == 'sp' and not explicit:
                code = None
            taxon = ('B. ' if abbreviated else 'Bacillus ') + match.group('species')
            snippet = text[max(0, match.start()-45):min(len(text), match.end()+65)].replace('\n', ' ')
            hits.append((taxon, code or '', snippet, source_field))
    return hits


def main():
    if not SOURCE.exists():
        raise SystemExit(f'Missing input: {SOURCE}')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    n_records = 0
    pmids_with_mentions = set()
    pmids_with_strain = set()
    counts = Counter()
    seen = set()
    with SOURCE.open(encoding='utf-8', newline='') as inp, OUT.open('w', encoding='utf-8', newline='') as out:
        reader = csv.DictReader(inp, delimiter='\t')
        required = {'pmid', 'doi', 'year', 'title', 'abstract', 'search_groups'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f'Missing input columns: {sorted(missing)}')
        writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        for rec in reader:
            n_records += 1
            for field in ('title', 'abstract'):
                for taxon, code, snippet, fieldname in mentions(rec.get(field) or '', field):
                    key = (rec['pmid'], taxon.lower(), code.lower(), fieldname)
                    if key in seen:
                        continue
                    seen.add(key)
                    pmids_with_mentions.add(rec['pmid'])
                    if code:
                        pmids_with_strain.add(rec['pmid'])
                    level = 'candidate_strain' if code else 'species_or_genus'
                    counts[level] += 1
                    writer.writerow({
                        'pmid': rec['pmid'], 'doi': rec['doi'], 'year': rec['year'],
                        'title': rec['title'], 'search_groups': rec['search_groups'],
                        'mention_text': snippet, 'species_mention': taxon,
                        'strain_candidate': code, 'mention_level': level,
                        'source_field': fieldname, 'verification_status': 'unverified_text_mention'
                    })
    report = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source': str(SOURCE), 'records_processed': n_records,
        'records_with_bacillus_mentions': len(pmids_with_mentions),
        'records_with_candidate_strain_mentions': len(pmids_with_strain),
        'mention_rows': sum(counts.values()), 'mention_types': dict(counts),
        'output': str(OUT),
        'limitations': 'Regex-based text mentions; abbreviated species names can be ambiguous; strain codes may be false positives. No eligibility or genome correspondence established.'
    }
    SUMMARY.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'Processed {n_records} records; wrote {sum(counts.values())} mentions.')
    print(f'Mentions: {OUT}\nSummary: {SUMMARY}')

if __name__ == '__main__':
    main()
