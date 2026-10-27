#!/usr/bin/env python3
"""Conservative triage of numeric mentions; never validate strain-level phenotypes.

Inputs: module 10 numeric mentions and module 11 publication-level candidates.
No numerical normalization, no imputation, no strain-to-value assignment.
"""
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NUM = ROOT / 'data/interim/pubmed_quantitative_mentions.tsv'
LINK = ROOT / 'data/interim/pubmed_strain_phenotype_candidates.tsv'
OUT = ROOT / 'data/interim/pubmed_numeric_endpoint_triage.tsv'
LINKOUT = ROOT / 'data/interim/pubmed_numeric_genome_context.tsv'
REPORT = ROOT / 'results/reports/pubmed_numeric_endpoint_triage_summary.json'

ENDPOINT = re.compile(r'\b(biosorp\w*|adsorp\w*|sorp\w*|remov\w*|uptake|bind\w*|biofilm\w*|biomass|exopolysaccharid\w*|extracellular polymeric|EPS|crystal violet|OD\s*\d{3}|matrix|production|yield|capacity|efficien\w*)\b', re.I)
CONDITIONS = re.compile(r'\b(initial|starting|supplied|added|treated|treatment|exposed|exposure|incubat\w*|cultiv\w*|medium|media|concentration|dose|dosage|temperature|shak\w*|contact time|after|before|at\s+\d+|for\s+\d+)\b', re.I)
OUTCOMES = re.compile(r'\b(remov\w*|biosorp\w*|adsorp\w*|sorp\w*|accumulat\w*|bind\w*|biofilm\w*|exopolysaccharid\w*|extracellular polymeric|EPS|yield|capacity|efficien\w*|production|increas\w*|decreas\w*|reduc\w*|enhanc\w*)\b', re.I)


def read(path, required):
    if not path.is_file():
        raise SystemExit(f'Missing input: {path}')
    with path.open(encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f'Missing columns in {path}: {sorted(missing)}')
        return list(reader)


def classify(r):
    unit = re.sub(r'\s+', '', r['unit_text'].lower()).replace('−','-')
    context = r['context_text'] or ''
    domains = set((r['domains'] or '').split(';'))
    # Candidate only; full sentence may contain several unrelated quantities.
    if unit in ('h','hour','hours','day','days','min','nm','mm','µm','μm'):
        return 'experimental_condition_or_method', 'time_wavelength_or_concentration_unit'
    if unit in ('mg/l','mgl-1','g/l','gl-1','µg/ml','μg/ml'):
        return 'ambiguous_numeric_context', 'solution_concentration_could_be_input_or_residual'
    if unit.startswith('od'):
        return 'ambiguous_numeric_context', 'optical_density_wavelength_or_measurement_unclear'
    if unit in ('mg/g','mgg-1'):
        if 'metal_sorption' in domains and OUTCOMES.search(context):
            return 'candidate_metal_sorption_endpoint', 'mass_per_mass_with_sorption_context'
        return 'ambiguous_numeric_context', 'mass_per_mass_without_sorption_endpoint_context'
    if unit == '%':
        if CONDITIONS.search(context) and not OUTCOMES.search(context):
            return 'experimental_condition_or_method', 'percentage_with_condition_context'
        if 'metal_sorption' in domains and OUTCOMES.search(context):
            return 'candidate_metal_removal_endpoint', 'percentage_with_sorption_context'
        if ('biofilm' in domains or 'eps_matrix' in domains) and OUTCOMES.search(context):
            return 'candidate_biofilm_eps_endpoint', 'percentage_with_biofilm_or_eps_context'
        return 'ambiguous_numeric_context', 'percentage_without_specific_endpoint_context'
    return 'ambiguous_numeric_context', 'unit_not_diagnostic'


def main():
    nums = read(NUM, ('pmid','doi','source_field','sentence_index','value_text','unit_text','context_text','domains','metals'))
    links = read(LINK, ('pmid','source_field','sentence_index','species_mention','strain_candidate','assembly_accession','biosample_accession'))
    by_sentence = defaultdict(list)
    for r in links:
        by_sentence[(r['pmid'],r['source_field'],r['sentence_index'])].append(r)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    fields = ['numeric_mention_id','pmid','doi','source_field','sentence_index','value_text','unit_text','domains','metals','triage_class','triage_reason','context_text','verification_status']
    linked_fields = fields + ['species_mention','strain_candidate','assembly_accession','biosample_accession','association_level']
    counts = Counter()
    linked_pmids = set()
    with OUT.open('w',encoding='utf-8',newline='') as fo, LINKOUT.open('w',encoding='utf-8',newline='') as fl:
        w = csv.DictWriter(fo,fieldnames=fields,delimiter='\t'); w.writeheader()
        wl = csv.DictWriter(fl,fieldnames=linked_fields,delimiter='\t'); wl.writeheader()
        for i, r in enumerate(nums,1):
            cls, reason = classify(r)
            row = {'numeric_mention_id':f'N{i:06d}',**{k:r.get(k,'') for k in ('pmid','doi','source_field','sentence_index','value_text','unit_text','domains','metals','context_text')},'triage_class':cls,'triage_reason':reason,'verification_status':'automated_unverified_not_strain_assigned'}
            w.writerow(row)
            counts[cls] += 1
            key = (r['pmid'],r['source_field'],r['sentence_index'])
            # Avoid multiplying numeric mentions by repeated evidence rows with identical genome/strain context.
            seen = set()
            for g in by_sentence.get(key,[]):
                ident = (g['species_mention'],g['strain_candidate'],g['assembly_accession'],g['biosample_accession'])
                if ident in seen: continue
                seen.add(ident)
                wl.writerow({**row,**{k:g[k] for k in ('species_mention','strain_candidate','assembly_accession','biosample_accession')},'association_level':'same_publication_only_not_strain_attributed'})
                counts['numeric_genome_context_rows'] += 1
                linked_pmids.add(r['pmid'])
    report = {'generated_utc':datetime.now(timezone.utc).isoformat(),'input_numeric_rows':len(nums),'input_publication_context_rows':len(links),'triage_class_rows':{k:counts[k] for k in sorted(counts) if k not in ('numeric_genome_context_rows',)},'numeric_genome_context_rows':counts['numeric_genome_context_rows'],'publications_with_numeric_and_strain_context':len(linked_pmids),'output':str(OUT.relative_to(ROOT)),'context_output':str(LINKOUT.relative_to(ROOT)),'limitations':'Heuristic sentence-level triage; numerical mentions are not validated experimental endpoints. Same-publication strain and genome context does not establish strain-to-measurement attribution. Do not pool units or treat candidate classes as verified labels.'}
    REPORT.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(f"Triaged {len(nums)} numeric mentions; linked context rows: {counts['numeric_genome_context_rows']}")
    print(f'Report: {REPORT.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
