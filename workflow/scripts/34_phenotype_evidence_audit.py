#!/usr/bin/env python3
"""Conservative, read-only phenotype evidence audit for bmetal.

Uses existing project TSVs; does not perform remote lookups or claim manual validation.
"""
import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'No TSV header: {path}')
        return list(reader), set(reader.fieldnames)


def write_tsv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def clean(value):
    return (value or '').strip()


def norm(value):
    return re.sub(r'\s+', ' ', clean(value)).casefold()


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, default=Path('.'))
    args = ap.parse_args()
    root = args.root.resolve()
    inp = {
        'integrated': root / 'data/processed/genome_phenotype_integrated.tsv',
        'verified': root / 'data/processed/verified_genome_phenotype_pairs.tsv',
        'strain': root / 'data/interim/pubmed_strain_phenotype_candidates.tsv',
        'qmax': root / 'data/processed/biosorption_qmax_observations_candidate.tsv',
    }
    for label, path in inp.items():
        if not path.is_file():
            raise FileNotFoundError(f'{label}: {path}')
    data = {label: read_tsv(path)[0] for label, path in inp.items()}
    fields = {label: read_tsv(path)[1] for label, path in inp.items()}
    required = {
        'integrated': {'pmid','species','strain','assembly_accession','link_status','evidence_text'},
        'verified': {'species','strain','assembly_accession','genome_identity_status','phenotype_status','model_eligibility'},
        'strain': {'pmid','species_mention','strain_candidate','assembly_accession','domains','association_level','evidence_text'},
        'qmax': {'pmid','species_strain_candidates','capacity_kind','evidence_text'},
    }
    for label, cols in required.items():
        missing = cols - fields[label]
        if missing:
            raise ValueError(f'{label}: missing columns {sorted(missing)}')

    # Each row is a *publication/strain/assembly candidate*, not a verified pair.
    candidates = {}
    for r in data['strain']:
        species, strain, assembly, pmid = (clean(r.get(k)) for k in ('species_mention','strain_candidate','assembly_accession','pmid'))
        if not species or not strain or not assembly or not pmid:
            continue
        domains = {x.strip() for x in clean(r.get('domains')).split(';') if x.strip()}
        if not domains.intersection({'biofilm','eps_matrix','metal_sorption'}):
            continue
        key = (norm(species), norm(strain), assembly, pmid)
        if key not in candidates:
            candidates[key] = {'species':species,'strain':strain,'assembly_accession':assembly,
                'pmid':pmid,'biosample_accession':'','domains':set(), 'evidence':[],
                'n_source_mentions':0,'numeric_mentions':set(), 'sources':set(),
                'source_link_statuses':set()}
        c = candidates[key]
        c['domains'].update(domains)
        c['n_source_mentions'] += 1
        c['sources'].add('pubmed_strain_phenotype_candidates')
        if clean(r.get('biosample_accession')):
            c['biosample_accession'] = clean(r['biosample_accession'])
        c['source_link_statuses'].add(clean(r.get('association_level')))
        mention = clean(r.get('numeric_mentions_in_sentence'))
        if mention:
            c['numeric_mentions'].add(mention)
        ev = clean(r.get('evidence_text'))
        if ev and ev not in c['evidence']:
            c['evidence'].append(ev)

    # Attach integrated observations only when the PMID, strain, species and accession
    # all match; no inference from publication-only or species-only matches.
    integrated_by_key = defaultdict(list)
    for r in data['integrated']:
        species, strain, assembly, pmid = (clean(r.get(k)) for k in ('species','strain','assembly_accession','pmid'))
        if all((species,strain,assembly,pmid)):
            integrated_by_key[(norm(species),norm(strain),assembly,pmid)].append(r)

    verified_by_key = defaultdict(list)
    for r in data['verified']:
        key = (norm(r.get('species')),norm(r.get('strain')),clean(r.get('assembly_accession')))
        verified_by_key[key].append(r)

    outrows = []
    for key, c in sorted(candidates.items(), key=lambda kv: kv[0]):
        obs = integrated_by_key.get(key, [])
        verified = verified_by_key.get(key[:3], [])
        statuses = {clean(r.get('link_status')) for r in obs if clean(r.get('link_status'))}
        c['source_link_statuses'].update(statuses)
        phenotype_evidence = [clean(r.get('evidence_text')) for r in obs if clean(r.get('evidence_text'))]
        all_ev = list(dict.fromkeys(c['evidence'] + phenotype_evidence))
        identity = any(clean(r.get('genome_identity_status')) == 'confirmed_exact_strain' for r in verified)
        verified_record = bool(verified) and identity
        # A registered exact strain match is not automatically a verified phenotype.
        decision = 'manual_primary_source_review'
        reason = 'publication_level_association_only; phenotype_not_verified'
        if verified_record:
            decision = 'registered_exact_strain_pending_phenotype_audit'
            reason = 'exact_strain_registered; phenotype_measurement_requires_primary_source_audit'
        elif obs:
            reason = 'integrated_candidate; identity_and_phenotype_require_verification'
        if len(c['domains']) == 1 and 'eps_matrix' in c['domains']:
            reason += '; eps_endpoint_may_not_measure_biofilm_or_biosorption'
        outrows.append({
            'species':c['species'], 'strain':c['strain'], 'assembly_accession':c['assembly_accession'],
            'biosample_accession':c['biosample_accession'], 'pmid':c['pmid'],
            'domains':';'.join(sorted(c['domains'])),
            'n_source_mentions':c['n_source_mentions'],
            'numeric_mentions': ' | '.join(sorted(c['numeric_mentions'])),
            'integrated_observation_rows':len(obs),
            'integrated_capacity_kinds':';'.join(sorted({clean(r.get('capacity_kind')) for r in obs if clean(r.get('capacity_kind'))})),
            'source_link_statuses':';'.join(sorted(x for x in c['source_link_statuses'] if x)),
            'registered_exact_strain': 'yes' if verified_record else 'no',
            'audit_decision':decision, 'audit_reason':reason,
            'evidence_excerpt': ' | '.join(all_ev[:2])[:1500],
            'source_tables':';'.join(sorted(c['sources'] | ({'genome_phenotype_integrated'} if obs else set()))),
        })

    outdir = root / 'results/functional/phenotype_audit'
    outdir.mkdir(parents=True, exist_ok=True)
    auditpath = outdir / 'phenotype_audit.tsv'
    fields_out = ['species','strain','assembly_accession','biosample_accession','pmid','domains',
        'n_source_mentions','numeric_mentions','integrated_observation_rows',
        'integrated_capacity_kinds','source_link_statuses','registered_exact_strain',
        'audit_decision','audit_reason','evidence_excerpt','source_tables']
    write_tsv(auditpath, outrows, fields_out)

    # Preserve original verified registry unchanged. Do not fabricate verified phenotype rows.
    # Export an exact copy with an explicit warning that phenotype audit is still pending.
    verifiedpath = outdir / 'phenotype_verified_final.tsv'
    verified_fields = list(read_tsv(inp['verified'])[1])
    # Keep input column order rather than unordered set.
    with inp['verified'].open(encoding='utf-8-sig', newline='') as f:
        verified_fields = csv.DictReader(f, delimiter='\t').fieldnames
    write_tsv(verifiedpath, data['verified'], verified_fields)

    summary = {
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'automated triage of existing evidence; NOT manual primary-source verification',
        'input_rows': {k:len(v) for k,v in data.items()},
        'candidate_publication_strain_assembly_rows':len(outrows),
        'distinct_candidate_species_strains':len({(norm(r['species']),norm(r['strain'])) for r in outrows}),
        'distinct_candidate_assemblies':len({r['assembly_accession'] for r in outrows}),
        'distinct_candidate_publications':len({r['pmid'] for r in outrows}),
        'audit_decisions':dict(Counter(r['audit_decision'] for r in outrows)),
        'registered_exact_strain_rows_in_verified_registry':len(data['verified']),
        'newly_verified_phenotypes_by_this_script':0,
        'verified_export_note':'Exact copy of prior registry, not newly phenotype-verified; consult phenotype_status and model_eligibility.',
        'limitations':['Publication-level co-occurrence is not phenotype assignment.',
                       'Automated evidence snippets cannot establish experimental comparability.',
                       'Strain aliases and generic strain names are not automatically merged.',
                       'Different metals, biomass states and units cannot be pooled without source review.'],
        'inputs_sha256':{k:sha256(v) for k,v in inp.items()},
        'outputs_sha256':{auditpath.name:sha256(auditpath),verifiedpath.name:sha256(verifiedpath)},
    }
    summarypath = outdir / 'phenotype_audit_summary.json'
    summarypath.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('inputs_sha256','outputs_sha256','limitations')},indent=2,ensure_ascii=False))
    print(f'Output directory: {outdir}')

if __name__ == '__main__':
    main()
