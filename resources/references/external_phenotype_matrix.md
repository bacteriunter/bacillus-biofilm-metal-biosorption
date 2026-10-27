# External phenotype-linked genome matrix

## Scope

Two genome assemblies were analyzed separately from the
638-genome main cohort:

- GCF_001858395.1 — AB211
- GCA_001880305.1 — AR156

Their genomic marker detections are incorporated into
manuscript Table 4 by workflow module 43.

## Input genomes

The original genome paths are recorded in:

data/processed/external_phenotype_fastas.txt

The genome FASTA files are external inputs and are not
distributed as part of the repository.

## Module 27: Homology searches

Software recorded:
MMseqs2 15-6f452+ds-2

Parameters:

- Search type: 2
- Sensitivity: 5.7
- E-value: 1e-5
- Threads: 2
- Reference proteins: 99

Command:

```bash
python workflow/scripts/27_search_external_phenotype_genomes.py \
  --fasta-list data/processed/external_phenotype_fastas.txt \
  --references data/reference/functional_panel/functional_reference_panel_merged.faa \
  --outdir results/functional/external_phenotypes/module27 \
  --threads 2 \
  --sensitivity 5.7 \
  --evalue 1e-5
```

Original execution metadata:

results/functional/external_phenotypes/module27/run_config.json

## Module 28: Functional marker matrix

Parameters:

- Maximum E-value: 1e-5
- Minimum query coverage: 0.7
- Minimum percent identity: 0

Command:

```bash
python workflow/scripts/28_build_functional_matrix.py \
  --alignments results/functional/external_phenotypes/module27 \
  --catalog data/reference/functional_panel/functional_reference_panel_merged.tsv \
  --reference-fasta data/reference/functional_panel/functional_reference_panel_merged.faa \
  --cohort data/processed/external_phenotype_fastas.txt \
  --outdir results/functional/external_phenotypes/module28 \
  --max-evalue 1e-5 \
  --min-qcov 0.7 \
  --min-pident 0
```

Original processing summary:

- Genomes: 2
- Reference proteins: 99
- Functional marker labels: 53
- Raw annotated hits: 966
- Filtered annotated hits: 718
- Status: PASS

The original output checksums are recorded in:

results/functional/external_phenotypes/module28/summary.json

## Manuscript integration

Module 43 combines this external matrix with the
638-genome main matrix to produce the phenotype-genome
evidence table.

The manuscript identifies this output as Table 4,
although the original output filename retains `table3`.

## Interpretation and limitations

Marker detection represents candidate sequence homology,
not experimentally confirmed gene function or orthology.

The external genomes were not incorporated into the
638-genome comparative cohort.

Reconstruction requires independently obtaining the
original genome assemblies and functional references.

The commands above document the recorded processing
configuration; independent reconstruction has not yet
been validated.
