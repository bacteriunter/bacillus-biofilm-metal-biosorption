#!/usr/bin/env python3
"""Select a reproducible Bacillus genome cohort from NCBI Datasets metadata.

No genome FASTA download, no implicit cap, no phenotype-driven exclusions.
Requires NCBI datasets CLI and Python >=3.10. Run --help for details.
"""
import argparse, csv, json, re, shutil, subprocess, sys, hashlib
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = ['assembly_accession','species','strain','biosample','bioproject','assembly_level','refseq_category','checkm_completeness','checkm_contamination','total_sequence_length','contig_n50','release_date','phenotype_candidate','decision','reason','selection_rank']

def get(obj,*keys):
    for k in keys:
        if not isinstance(obj,dict): return None
        obj=obj.get(k)
    return obj

def first(*values):
    return next((x for x in values if x is not None and x != ''),None)

def num(x):
    try:return float(x)
    except (ValueError,TypeError):return None

def clean(x):return str(x if x is not None else '').strip()

def species_of(name):
    # Retain explicit Bacillus binomials; preserve bracketed names
    # as distinct taxonomic labels rather than silently renaming them.
    parts = name.strip().split()
    if len(parts) < 2 or parts[0] != 'Bacillus':
        return ''
    epithet = parts[1]
    if epithet.lower() in ('sp.', 'sp', 'spp.', 'spp'):
        return ''
    if re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*', epithet):
        return 'Bacillus ' + epithet
    if re.fullmatch(r'\\[[A-Za-z][A-Za-z0-9_-]*\\]', epithet):
        return 'Bacillus ' + epithet
    return ''

def load_phenotypes(path):
    identifiers=set(); rows=[]
    if not path.exists():return identifiers,rows
    with path.open(encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            rows.append(r)
            for k in ('assembly_accession','verified_assembly','assembly_candidates','biosample_accession','verified_biosample','biosample_candidates'):
                for a in re.findall(r'\b(?:GC[AF]_\d+\.\d+|SAM[NED]\d+)\b',r.get(k,'') or ''):identifiers.add(a)
    return identifiers,rows

def parse_record(r,ids):
    org=r.get('organism') or {}; info=r.get('assembly_info') or {}; stats=r.get('assembly_stats') or {}; check=info.get('checkm_info') or r.get('checkm_info') or {}
    name=clean(org.get('organism_name'))
    strain=first(get(org,'infraspecific_names','strain'),get(info,'biosample','strain'))
    biosample=first(get(info,'biosample','accession'),r.get('biosample_accession'))
    accession=clean(r.get('accession'))
    return dict(assembly_accession=accession,species=species_of(name),strain=clean(strain),biosample=clean(biosample),bioproject=clean(info.get('bioproject_accession')),assembly_level=clean(info.get('assembly_level')),refseq_category=clean(info.get('refseq_category')),checkm_completeness=clean(first(check.get('completeness'),check.get('checkm_completeness'))),checkm_contamination=clean(first(check.get('contamination'),check.get('checkm_contamination'))),total_sequence_length=clean(first(stats.get('total_sequence_length'),stats.get('total_sequence_len'))),contig_n50=clean(first(stats.get('contig_n50'),stats.get('contig_n50_length'))),release_date=clean(info.get('release_date')),phenotype_candidate='yes' if accession in ids or biosample in ids else 'no',decision='',reason='',selection_rank='')

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project',default='.')
    ap.add_argument('--metadata',help='Existing NCBI Datasets JSONL; if absent, query NCBI for all RefSeq Bacillus genomes')
    ap.add_argument('--max-per-species',type=int,default=20,help='Proposed species cap; selection by deterministic quality ranking (default 20)')
    ap.add_argument('--min-completeness',type=float,default=95)
    ap.add_argument('--max-contamination',type=float,default=5)
    ap.add_argument('--allow-missing-checkm',action='store_true',help='Include unassessed assemblies as provisional, clearly labeled')
    args=ap.parse_args()
    if args.max_per_species<1:ap.error('--max-per-species must be positive')
    root=Path(args.project); raw=root/'data/raw/ncbi_bacillus_refseq_metadata.jsonl'
    if args.metadata:raw=Path(args.metadata)
    else:
        if not shutil.which('datasets'):sys.exit('ERROR: NCBI datasets CLI missing')
        raw.parent.mkdir(parents=True,exist_ok=True)
        print('Querying all RefSeq Bacillus assemblies. This can take time; no record cap is applied.',flush=True)
        cmd=['datasets','summary','genome','taxon','1386','--assembly-source','RefSeq','--as-json-lines']
        with raw.open('w',encoding='utf-8') as out:
            result=subprocess.run(cmd,stdout=out)
        if result.returncode:sys.exit(f'NCBI query failed ({result.returncode}); incomplete file {raw} must not be used')
    ids, phenotypes=load_phenotypes(root/'data/processed/genome_phenotype_integrated.tsv')
    for r in phenotypes:
        for k in ('assembly_accession','verified_assembly','biosample_accession'):
            if r.get(k):ids.add(r[k].strip())
    records={}; malformed=0
    with raw.open(encoding='utf-8') as f:
        for line in f:
            if not line.strip():continue
            try:
                data=json.loads(line)
                if 'reports' in data:
                    items=data['reports']
                else:items=[data]
                for item in items:
                    if not isinstance(item,dict) or not item.get('accession'):continue
                    rec=parse_record(item,ids)
                    if rec['assembly_accession'].startswith('GCF_'):records[rec['assembly_accession']]=rec
            except (ValueError,TypeError):malformed+=1
    if not records:sys.exit('ERROR: no GCF assemblies parsed; check JSONL format and query')
    by_species=defaultdict(list); rows=[]
    for rec in records.values():
        c=num(rec['checkm_completeness']); t=num(rec['checkm_contamination'])
        if not rec['species']:
            rec.update(decision='excluded',reason='species_not_resolved')
        elif c is None or t is None:
            if args.allow_missing_checkm:rec.update(decision='candidate',reason='quality_unassessed')
            else:rec.update(decision='pending_quality',reason='checkm_metrics_missing')
        elif c<args.min_completeness or t>args.max_contamination:
            rec.update(decision='excluded',reason='quality_threshold_failed')
        else:rec.update(decision='candidate',reason='quality_pass')
        if rec['decision']=='candidate':by_species[rec['species']].append(rec)
        rows.append(rec)
    # Preserve all exact-identifier phenotype candidates in a separate priority stratum;
    # do not assert that their phenotype link is verified.
    selected=[]
    for sp,items in sorted(by_species.items()):
        def score(r):
            level={'Complete Genome':4,'Chromosome':3,'Scaffold':2,'Contig':1}.get(r['assembly_level'],0)
            return (r['phenotype_candidate']!='yes',-level,-(num(r['checkm_completeness']) or 0),num(r['checkm_contamination']) or 0,-(num(r['contig_n50']) or 0),r['assembly_accession'])
        items.sort(key=score)
        priority=[r for r in items if r['phenotype_candidate']=='yes']
        other=[r for r in items if r['phenotype_candidate']!='yes']
        chosen=priority+other[:max(0,args.max_per_species-len(priority))]
        for rank,r in enumerate(chosen,1):
            r['decision']='selected';r['reason']='phenotype_identifier_priority' if r['phenotype_candidate']=='yes' else 'species_quality_rank';r['selection_rank']=str(rank)
        for r in items:
            if r['decision']=='candidate':r.update(decision='not_selected',reason='species_cap')
        selected.extend(chosen)
    rows.sort(key=lambda r:(r['species'],r['assembly_accession']))
    selected.sort(key=lambda r:(r['species'],int(r['selection_rank'])))
    processed=root/'data/processed';reports=root/'results/reports'
    write_tsv(processed/'genome_sampling_manifest.tsv',rows,FIELDS)
    write_tsv(processed/'genome_cohort_selected.tsv',selected,FIELDS)
    # Link phenotypes only by explicit exact Assembly or BioSample identifier; no PMID-only links.
    assembly_map={r['assembly_accession']:r for r in selected};bio_map=defaultdict(list)
    for r in selected:
        if r['biosample']:bio_map[r['biosample']].append(r)
    linked=[]
    for p in phenotypes:
        hits={}
        for key,value in p.items():
            if key in ('assembly_accession','verified_assembly','assembly_candidates','biosample_accession','biosample_candidates'):
                for a in re.findall(r'\b(?:GC[AF]_\d+\.\d+|SAM[NED]\d+)\b',value or ''):
                    if a in assembly_map:hits[a]=assembly_map[a]
                    for r in bio_map.get(a,[]):hits[r['assembly_accession']]=r
        for r in hits.values():
            linked.append(dict(p,selected_assembly=r['assembly_accession'],selected_biosample=r['biosample'],link_basis='matching_accession_candidate_not_verified'))
    columns=list(dict.fromkeys([k for p in phenotypes for k in p]+['selected_assembly','selected_biosample','link_basis']))
    write_tsv(processed/'phenotype_linked_subset.tsv',linked,columns)
    report={'metadata_file':str(raw),'metadata_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'parsed_refseq_assemblies':len(rows),'selected_assemblies':len(selected),'species_selected':len(by_species),'status_counts':dict(Counter(r['decision'] for r in rows)),'phenotype_identifier_priority_selected':sum(r['phenotype_candidate']=='yes' for r in selected),'phenotype_link_rows_candidate_only':len(linked),'malformed_json_lines':malformed,'parameters':vars(args),'limitations':['No ANI dereplication: this is a provisional quality-ranked cohort, not an ANI-dereplicated final cohort.','Missing CheckM metadata remain pending unless explicitly allowed.','Phenotype links are accession matches, not verified genotype-phenotype observations.','Species cap is a configurable design parameter, not a statistically optimized sample size.']}
    reports.mkdir(parents=True,exist_ok=True)
    (reports/'genome_cohort_selection_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
