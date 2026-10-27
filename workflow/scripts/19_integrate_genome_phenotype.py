#!/usr/bin/env python3
"""Integrate existing bmetal evidence without promoting candidate links to verified pairs.
Python >=3.10, standard library only. Run from project root.
"""
import argparse,csv,json,re,hashlib
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone

FIELDS=['record_id','pmid','doi','species','strain','metal','capacity_value','capacity_unit','capacity_kind','material_category','biomass_basis','source_file','evidence_text','assembly_accession','biosample_accession','genome_fasta','genome_identity_status','phenotype_status','link_status','model_eligibility','link_basis','source_tables']
INPUTS={
 'qmax':'data/processed/biosorption_qmax_observations_candidate.tsv',
 'fulltext':'data/processed/fulltext_biosorption_capacity_candidates.tsv',
 'links':'data/processed/genome_phenotype_publication_links.tsv',
 'strain':'data/interim/pubmed_strain_phenotype_candidates.tsv',
 'assembly':'data/interim/ncbi_assembly_candidates.tsv',
 'biosample':'data/interim/ncbi_biosample_resolution.tsv',
 'lookups':'data/processed/strain_assembly_lookup.tsv',
 'confirmed':'data/processed/verified_genome_phenotype_pairs.tsv',
}
def load(p):
 if not p.is_file(): return []
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def get(r,*keys):
 for k in keys:
  v=(r.get(k) or '').strip()
  if v:return v
 return ''
def norm(x):return re.sub(r'[^a-z0-9]','',x.lower())
def tokens(x):return [z.strip() for z in re.split(r'[;|]',x or '') if z.strip()]
def category(text):
 t=(text or '').lower()
 if re.search(r'\b(spore|spores|endospore|endospores|esporas?)\b',t):return 'spores'
 if re.search(r'\b(isolated eps|purified eps|extracted eps|extracellular polymeric substances|exopolysaccharides?|isolated exopolymers?|eps aislados?)\b',t):return 'isolated_eps'
 if re.search(r'\b(biosurfactant|bio-surfactant|microcapsule|alginate.chitosan|immobiliz|encapsulat)\b',t):return 'other_or_composite'
 if re.search(r'\b(live cells?|dead cells?|living cells?|non.living cells?|dry biomass|wet biomass|bacterial biomass|cell biomass|biomasa celular|whole.cells?)\b',t):return 'cell_biomass'
 return 'undetermined'
def choose_material(r):
 direct=get(r,'material_category','material_type','material_flag')
 if direct.lower() in ('spores','isolated_eps','cell_biomass','other_or_composite'):return direct.lower()
 return category(' '.join([get(r,'context','evidence_text'),direct]))
def assemble_row(r,origin):
 pmid=get(r,'pmid');doi=get(r,'doi')
 raw=get(r,'species_strain_candidates','species_mention','species','verified_strain')
 strain=get(r,'verified_strain','strain_candidate','strain')
 species=get(r,'species','species_mention')
 if not species and raw:
  m=re.search(r'\b(Bacillus\s+(?:[A-Za-z][\w-]*|sp\.))',raw,re.I)
  if m:species=m.group(1)
 if not strain and raw and species and norm(raw)!=norm(species):
  tail=raw[len(species):].strip(' ,;') if raw.lower().startswith(species.lower()) else ''
  if tail and not tail.lower().startswith(('strain ','subsp.')):strain=tail.split(';')[0].strip()
 val=get(r,'capacity_value','value','verified_capacity_mg_g_dry')
 kind=get(r,'capacity_kind','triage_class') or ('unclassified_fulltext_mention' if origin=='fulltext' else 'unknown')
 unit=get(r,'capacity_unit','unit_original') or ('mg/g dry' if get(r,'verified_capacity_mg_g_dry') else '')
 basis=get(r,'biomass_basis') or 'unspecified'
 evidence=get(r,'evidence_text','context')
 metal=get(r,'metal','metal_mention')
 assembly=get(r,'verified_assembly','assembly_accession','assembly_candidates')
 biosample=get(r,'biosample_accession','biosample_candidates')
 material=choose_material(r)
 rid=hashlib.sha256('|'.join([origin,pmid,doi,raw,val,unit,kind,evidence,get(r,'sentence_index')]).encode()).hexdigest()[:16]
 return dict(record_id=rid,pmid=pmid,doi=doi,species=species,strain=strain,metal=metal,capacity_value=val,capacity_unit=unit,capacity_kind=kind,material_category=material,biomass_basis=basis,source_file=get(r,'source_file'),evidence_text=evidence,assembly_accession=assembly,biosample_accession=biosample,genome_fasta=get(r,'genome_fasta'),genome_identity_status=get(r,'genome_identity_status') or 'unverified',phenotype_status=get(r,'phenotype_status','evidence_status') or 'candidate',link_status='unlinked',model_eligibility='not_assessed',link_basis='',source_tables=origin)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',default='.');a=ap.parse_args();root=Path(a.root)
 data={k:load(root/v) for k,v in INPUTS.items()}
 if not data['qmax'] and not data['fulltext'] and not data['confirmed']:raise SystemExit('ERROR: no phenotype evidence input files found')
 rows=[];seen=set()
 # qmax is consolidated primary source; fulltext only fills publications not in qmax to avoid duplicate extraction
 for r in data['qmax']:
  x=assemble_row(r,'qmax');rows.append(x);seen.add(x['pmid'])
 for r in data['fulltext']:
  if get(r,'pmid') in seen:continue
  rows.append(assemble_row(r,'fulltext'))
 # Explicitly registered records are preserved, never silently merged with extracted candidates.
 for r in data['confirmed']:
  x=assemble_row(r,'registered');x['material_category']=choose_material(r)
  x['capacity_kind']='registered_reported_value';x['link_status']='registered_exact_strain_genome' if get(r,'genome_identity_status')=='confirmed_exact_strain' else 'registered_unverified'
  x['link_basis']='user_registered_record';x['model_eligibility']=get(r,'model_eligibility') or 'not_assessed'
  rows.append(x)
 # Candidate indices: match only same PMID + explicit strain; never publication alone.
 idx=defaultdict(list)
 for source in ('links','strain'):
  for r in data[source]:
   pmid=get(r,'pmid');s=get(r,'strain_candidate','strain','verified_strain');acc=get(r,'assembly_accession','assembly_candidates');bio=get(r,'biosample_accession','biosample_candidates')
   if pmid and s and (acc or bio):idx[(pmid,norm(s))].append((acc,bio,source))
 # Assembly metadata from lookup is only used if taxon is Bacillus AND exact strain name matches.
 lookup=defaultdict(list)
 for r in data['lookups']:
  organism=get(r,'ncbi_organism')
  if not re.search(r'\bBacillus\b',organism,re.I):continue
  if get(r,'match_status')!='exact_identifier_in_metadata':continue
  s=get(r,'requested_strain');acc=get(r,'assembly_accession');bio=get(r,'biosample_accession')
  if s and acc:lookup[norm(s)].append((acc,bio,'lookup_exact_identifier_candidate'))
 for x in rows:
  if x['source_tables']=='registered':continue
  if x['assembly_accession'] or x['biosample_accession']:
   x['link_status']='candidate_in_source';x['link_basis']='extracted_source_candidate';continue
  if not x['strain']:continue
  candidates=idx.get((x['pmid'],norm(x['strain'])),[])
  if not candidates:candidates=lookup.get(norm(x['strain']),[])
  uniq=sorted(set((aa,bb) for aa,bb,_ in candidates))
  if len(uniq)==1:
   x['assembly_accession'],x['biosample_accession']=uniq[0]
   x['link_status']='candidate_exact_strain_name';x['link_basis']=';'.join(sorted(set(c for _,_,c in candidates)))
  elif len(uniq)>1:
   x['link_status']='ambiguous_multiple_genomes';x['link_basis']=';'.join(sorted(set(c for _,_,c in candidates)))
 # No candidate is promoted to verified, even when its identifier matches.
 out=root/'data/processed/genome_phenotype_integrated.tsv';out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader();w.writerows(rows)
 report={'generated_utc':datetime.now(timezone.utc).isoformat(),'inputs':{k:{'path':v,'rows':len(data[k]),'exists':(root/v).is_file()} for k,v in INPUTS.items()},'output':str(out.relative_to(root)),'integrated_rows':len(rows),'material_counts':dict(Counter(r['material_category'] for r in rows)),'link_status_counts':dict(Counter(r['link_status'] for r in rows)),'capacity_kind_counts':dict(Counter(r['capacity_kind'] for r in rows)),'note':'Candidates are not verified genome-phenotype pairs. EPS, spores, cellular biomass, and other materials remain separate. No model eligibility is inferred.'}
 rep=root/'results/reports/genome_phenotype_integration_summary.json';rep.parent.mkdir(parents=True,exist_ok=True);rep.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
 print(f'Integrated {len(rows)} records -> {out}');print(f'Report -> {rep}')
 print('Material:',report['material_counts']);print('Links:',report['link_status_counts'])
if __name__=='__main__':main()
