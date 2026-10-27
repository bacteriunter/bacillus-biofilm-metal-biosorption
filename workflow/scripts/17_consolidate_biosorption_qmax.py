#!/usr/bin/env python3
"""Consolidate full-text biosorption capacity candidates without inventing phenotype-genome links.

Input: data/processed/fulltext_biosorption_capacity_candidates.tsv (module 16)
Output: candidate observations and an evidence-preserving report.
Python standard library only.
"""
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

SRC = Path('data/processed/fulltext_biosorption_capacity_candidates.tsv')
OUT = Path('data/processed/biosorption_qmax_observations_candidate.tsv')
REPORT = Path('results/reports/biosorption_qmax_consolidation_summary.json')

FIELDS = ['pmid','doi','species_strain_candidates','metal','capacity_mg_g','capacity_kind',
          'biomass_state','biomass_basis','material_flag','assembly_candidates',
          'biosample_candidates','source_file','evidence_sentence_indices','evidence_text',
          'supporting_mentions','screening_status','verification_notes']


def norm(s):
    return re.sub(r'\s+', ' ', (s or '').strip())


def biomass_state(text, value):
    t = text.lower()
    # Resolve specifically paired Qmax values, e.g. '12.94 and 20.35 mg/g for live and dead cells respectively'.
    pattern = re.search(r'(\d+(?:\.\d+)?)\s*(?:and|&)\s*(\d+(?:\.\d+)?)\s*mg\s*/\s*g\s+for\s+(live|dead|living|nonliving|non-living|inactivated)\s+and\s+(live|dead|living|nonliving|non-living|inactivated)\s+(?:cells|biomass)', t)
    if pattern:
        first, second = float(pattern.group(1)), float(pattern.group(2))
        if abs(value-first)<1e-7:
            return 'live' if pattern.group(3) in ('live','living') else 'dead'
        if abs(value-second)<1e-7:
            return 'live' if pattern.group(4) in ('live','living') else 'dead'
    # If both treatments appear but cannot attribute value, leave unresolved.
    has_live = bool(re.search(r'\b(?:live|living)\s+(?:cells|biomass)\b', t))
    has_dead = bool(re.search(r'\b(?:dead|inactivated|heat.killed|non.living)\s+(?:cells|biomass)\b', t))
    if has_live and has_dead:
        return 'mixed_unattributed'
    if has_live:
        return 'live_context_only'
    if has_dead:
        return 'dead_context_only'
    return 'unspecified'


def capacity_kind(text):
    t = text.lower()
    if re.search(r'\b(?:q\s*max|qmax|maximum (?:biosorption|adsorption|sorption) capacity|langmuir (?:maximum|capacity))\b', t):
        return 'qmax_candidate'
    if re.search(r'\b(?:initial (?:concentration|ratio)|dosage|dose|added|applied|feed concentration)\b', t):
        return 'condition_or_dose_candidate'
    if re.search(r'\b(?:uptake|adsorption|biosorption|sorption)\b', t):
        return 'observed_uptake_candidate'
    return 'unclassified_mg_g'


def mass_basis(text, existing):
    t = text.lower()
    if re.search(r'\b(?:dry weight|dry biomass|dry cell weight|dry mass|dried biomass|lyophiliz(?:ed|ation))\b', t):
        return 'dry_context_only'
    if re.search(r'\b(?:fresh weight|wet weight|wet biomass|fresh biomass)\b', t):
        return 'wet_context_only'
    return norm(existing) if existing and existing != 'unspecified' else 'unconfirmed'


def main():
    if not SRC.exists():
        raise SystemExit(f'Missing input: {SRC}')
    with SRC.open(encoding='utf-8', newline='') as f:
        rd = csv.DictReader(f, delimiter='\t')
        required = {'pmid','doi','value','unit_original','context','sentence_index','species_strain_candidates'}
        missing = required - set(rd.fieldnames or [])
        if missing:
            raise SystemExit('Missing columns: ' + ', '.join(sorted(missing)))
        rows = list(rd)

    groups = defaultdict(list)
    reasons = Counter()
    for r in rows:
        unit = norm(r.get('unit_original','')).lower().replace(' ','')
        if unit not in ('mg/g','mgg-1','mgg(-1)','mgg⁻¹'):
            reasons['other_unit'] += 1
            continue
        try:
            value = float(r['value'])
        except (TypeError, ValueError):
            reasons['invalid_value'] += 1
            continue
        if value <= 0:
            reasons['nonpositive_value'] += 1
            continue
        ctx = norm(r.get('context'))
        kind = capacity_kind(ctx)
        state = biomass_state(ctx, value)
        basis = mass_basis(ctx, r.get('biomass_basis',''))
        # Never silently merge unknown and known biomass states.
        key = (norm(r['pmid']), norm(r.get('doi')), norm(r.get('species_strain_candidates')),
               norm(r.get('metal_mention')).lower(), value, kind, state, basis,
               norm(r.get('material_flag')), norm(r.get('assembly_candidates')),
               norm(r.get('biosample_candidates')))
        groups[key].append(r)
        reasons[kind] += 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    results = []
    for key, supports in sorted(groups.items()):
        (pmid,doi,strain,metal,value,kind,state,basis,material,assembly,biosample) = key
        contexts = list(dict.fromkeys(norm(r.get('context')) for r in supports))
        indices = list(dict.fromkeys(norm(r.get('sentence_index')) for r in supports))
        sources = list(dict.fromkeys(norm(r.get('source_file')) for r in supports))
        notes = []
        if kind != 'qmax_candidate':
            notes.append('Not explicitly identified as Qmax in local context')
        if basis != 'dry_context_only':
            notes.append('Dry biomass denominator not established')
        if state in ('unspecified','mixed_unattributed'):
            notes.append('Biomass state not attributed to this value')
        if not strain:
            notes.append('Strain not identified')
        if not assembly:
            notes.append('No assembly candidate in source table')
        notes.append('Publication-level strain/genome references are NOT verified measurement links')
        results.append(dict(zip(FIELDS, [pmid,doi,strain,metal,f'{value:g}',kind,state,basis,material,
                         assembly,biosample,'; '.join(sources),'; '.join(indices),' || '.join(contexts),
                         str(len(supports)),'requires_primary_source_verification','; '.join(notes)])))
    with OUT.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    summary = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'input_rows': len(rows), 'distinct_candidate_records': len(results),
        'capacity_kind_counts': dict(Counter(r['capacity_kind'] for r in results)),
        'qmax_candidate_rows': sum(r['capacity_kind']=='qmax_candidate' for r in results),
        'qmax_with_dry_context_rows': sum(r['capacity_kind']=='qmax_candidate' and r['biomass_basis']=='dry_context_only' for r in results),
        'qmax_with_assembly_candidate_rows': sum(r['capacity_kind']=='qmax_candidate' and bool(r['assembly_candidates']) for r in results),
        'unique_pmids': len({r['pmid'] for r in results}),
        'verified_genome_phenotype_pairs': 0,
        'output': str(OUT),
        'limitations': 'Conservative sentence-level consolidation only. Qmax, dry-mass basis, strain, metal, and genome must be checked against methods/tables; no model-ready observations are certified.'
    }
    REPORT.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(f'Consolidated {len(rows)} mentions into {len(results)} candidate records.')
    print(f'Qmax candidates: {summary["qmax_candidate_rows"]}; report: {REPORT}')

if __name__ == '__main__':
    main()
