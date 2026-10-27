#!/usr/bin/env python3
"""Consolidate MMseqs translated-search hits; no new alignment is run.
Outputs raw and threshold-filtered hits, genome x gene and genome x group matrices.
"""
import argparse,csv,hashlib,json,re,sys
from collections import defaultdict,Counter
from pathlib import Path

FIELDS=[
    'query','target','pident','alnlen',
    'qstart','qend','qlen','tstart','tend',
    'evalue','bits','qcov','tcov'
]

def rows(path):
    with path.open(newline='') as fh:
        return list(csv.DictReader(fh,delimiter='\t'))

def fasta(path):
    seq={}; key=None; parts=[]
    with path.open() as fh:
        for line in fh:
            line=line.strip()
            if line.startswith('>'):
                if key: seq[key]=''.join(parts)
                key=line[1:].split()[0].split('|')[0];parts=[]
            elif line: parts.append(line)
        if key: seq[key]=''.join(parts)
    return seq

def write(path,cols,data):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(data)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--alignments',default='results/functional/module27')
    p.add_argument('--catalog',default='data/reference/functional_panel/functional_reference_panel_merged.tsv')
    p.add_argument('--reference-fasta',default='data/reference/functional_panel/functional_reference_panel_merged.faa')
    p.add_argument('--cohort',default='data/processed/functional_cohort_638_fastas.txt')
    p.add_argument('--outdir',default='results/functional/module28')
    p.add_argument('--max-evalue',type=float,default=1e-5)
    p.add_argument('--min-qcov',type=float,default=0.7,help='Query coverage fraction 0..1; provisional operational filter')
    p.add_argument('--min-pident',type=float,default=0,help='Identity in percent; zero means no identity filter')
    a=p.parse_args()
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    if not 0<=a.min_qcov<=1: p.error('--min-qcov must be 0..1')
    refs=rows(Path(a.catalog)); refseq=fasta(Path(a.reference_fasta))
    # Multiple annotations may share an accession: preserve every distinct annotation.
    ann=defaultdict(list)
    for r in refs:
        acc=(r.get('accession') or r.get('uniprot_accession') or r.get('protein_accession') or '').split('|')[0]
        gene=r.get('gene') or r.get('query_gene') or ''
        group=r.get('group') or r.get('family') or ''
        if not acc or not gene or not group: raise ValueError('Catalog needs accession, gene and group; missing in '+str(r))
        if (group,gene) not in ann[acc]: ann[acc].append((group,gene))
    if not ann: raise ValueError('Empty reference catalog')
    missing=set(ann)-set(refseq)
    if missing: raise ValueError('Catalog accessions missing from FASTA: '+str(sorted(missing)[:10]))
    paths=[x.strip() for x in Path(a.cohort).read_text().splitlines() if x.strip()]
    genomes=[Path(x).stem for x in paths]
    if len(genomes)!=len(set(genomes)): raise ValueError('Duplicate genome accession in cohort list')
    if len(genomes)!=638: print(f'NOTE: cohort has {len(genomes)} genomes (not 638)',file=sys.stderr)
    alnroot=Path(a.alignments)
    hits=[]; accepted=[]; missing_files=[]; bad=[]; unmatched=set()
    for genome in genomes:
        matches=list(alnroot.rglob(genome+'.tsv'))
        # Some versions store hits with another suffix; fail instead of guessing.
        if len(matches)!=1:
            missing_files.append({'genome':genome,'matches':len(matches)});continue
        file=matches[0]
        with file.open() as f:
            for n,line in enumerate(f,1):
                if not line.strip() or line.startswith('#'):continue
                parts=line.rstrip('\n').split('\t')
                if parts[0]=='query':continue
                if len(parts)<13:
                    bad.append(f'{file}:{n} expected >=13 columns, got {len(parts)}');continue
                raw=dict(zip(FIELDS,parts[:13]))
                acc=raw['query'].split('|')[0].split()[0]
                if acc not in ann:
                    unmatched.add(acc);continue
                try:
                    ident=float(raw['pident']); alnlen=int(raw['alnlen']); qcov=float(raw['qcov']); tcov=float(raw['tcov']); ev=float(raw['evalue']); bits=float(raw['bits'])
                except ValueError:
                    bad.append(f'{file}:{n} nonnumeric alignment value');continue
                # MMseqs qcov is a fraction. If percent is present, reject rather than silently rescale.
                if not (0<=qcov<=1 and 0<=tcov<=1):
                    bad.append(f'{file}:{n} coverage outside 0..1');continue
                passed=ev<=a.max_evalue and qcov>=a.min_qcov and ident>=a.min_pident
                for group,gene in ann[acc]:
                    r={'genome':genome,'reference_accession':acc,'group':group,'gene':gene,'target':raw['target'],
                       'pident':ident,'alnlen':alnlen,'qcov':qcov,'tcov':tcov,'evalue':ev,'bits':bits,
                       'passes_filter':int(passed)}
                    hits.append(r)
                    if passed: accepted.append(r)
    if missing_files or bad or unmatched:
        errors={'missing_or_ambiguous_alignment_files':missing_files[:30], 'bad_rows':bad[:30], 'unknown_reference_accessions':sorted(unmatched)[:30]}
        (out/'input_errors.json').write_text(json.dumps(errors,indent=2)+'\n')
        raise RuntimeError('Input mismatch: inspect '+str(out/'input_errors.json'))
    hitcols=['genome','reference_accession','group','gene','target','pident','alnlen','qcov','tcov','evalue','bits','passes_filter']
    write(out/'all_hits.tsv',hitcols,hits)
    write(out/'filtered_hits.tsv',hitcols,accepted)
    keys=sorted({pair for pairs in ann.values() for pair in pairs})
    # Count distinct target identifiers per gene; this is not an ortholog copy-number estimate.
    targets=defaultdict(set)
    for r in accepted: targets[(r['genome'],r['group'],r['gene'])].add(r['target'])
    columns=['genome']+[g+'::'+gene for g,gene in keys]
    presence=[]; counts=[]
    for genome in genomes:
        presence.append({'genome':genome,**{g+'::'+gene:int(bool(targets[(genome,g,gene)])) for g,gene in keys}})
        counts.append({'genome':genome,**{g+'::'+gene:len(targets[(genome,g,gene)]) for g,gene in keys}})
    write(out/'gene_presence.tsv',columns,presence)
    write(out/'gene_target_counts.tsv',columns,counts)
    groups=sorted({g for g,_ in keys}); gp=[]
    for genome in genomes:
        gp.append({'genome':genome,**{g:int(any(targets[(genome,gg,gene)] for gg,gene in keys if gg==g)) for g in groups}})
    write(out/'group_presence.tsv',['genome']+groups,gp)
    summary={'status':'PASS','genomes':len(genomes),'reference_sequences':len(refseq),'annotated_accessions':len(ann),
             'gene_group_labels':len(keys),'raw_annotated_hits':len(hits),'filtered_annotated_hits':len(accepted),
             'genomes_with_filtered_hits':len({r['genome'] for r in accepted}),
             'filter':{'max_evalue':a.max_evalue,'min_query_coverage':a.min_qcov,'min_percent_identity':a.min_pident},
             'interpretation':'sequence homology candidates; target counts are not validated gene copy counts',
             'files':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in out.glob('*.tsv')}}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='files'},indent=2))

if __name__=='__main__':main()
