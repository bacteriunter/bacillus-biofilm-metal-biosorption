# Bacillus biofilm and metal biosorption project
# Reproducible manuscript figures and tables
#
# Scientific analyses upstream of these outputs
# are preserved in workflow/scripts/.
#
# Run from the project root.

PYTHON = "python"

rule all:
    input:
        "results/functional/module30/summary.json",
        "results/functional/module31/summary.json",
        "results/functional/module33/summary.json",
        "results/figures/module37/figure_marker_associations.pdf",
        "results/figures/module38/figure1_functional_distribution.pdf",
        "results/figures/module39/figure2_variance_partition.pdf",
        "results/figures/module40/figure4_phytoremediation_framework.pdf",
        "results/tables/module41/table1_genomic_dataset.tsv",
        "results/tables/module42/table2_functional_variability.tsv",
        "results/tables/module44/table4_marker_associations.tsv",
        "results/tables/module43/table3_phenotype_genome_evidence.tsv"

rule figure1:
    input:
        "results/functional/module29/species_gene_prevalence.tsv"
    output:
        "results/figures/module38/figure1_functional_distribution.pdf",
        "results/figures/module38/figure1_functional_distribution.png",
        "results/figures/module38/figure1_functional_distribution.svg"
    shell:
        "{PYTHON} workflow/scripts/38_plot_functional_distribution.py"

rule figure2:
    input:
        "results/functional/module32/gene_variance_partition.tsv"
    output:
        "results/figures/module39/figure2_variance_partition.pdf",
        "results/figures/module39/figure2_variance_partition.png",
        "results/figures/module39/figure2_variance_partition.svg"
    shell:
        "{PYTHON} workflow/scripts/39_plot_variance_partition.py"

rule figure3:
    input:
        "results/statistics/module36/marker_associations.tsv"
    output:
        "results/figures/module37/figure_marker_associations.pdf",
        "results/figures/module37/figure_marker_associations.png",
        "results/figures/module37/figure_marker_associations.svg"
    shell:
        "{PYTHON} workflow/scripts/37_plot_marker_associations.py"

rule figure4:
    output:
        "results/figures/module40/figure4_phytoremediation_framework.pdf",
        "results/figures/module40/figure4_phytoremediation_framework.png",
        "results/figures/module40/figure4_phytoremediation_framework.svg"
    shell:
        "{PYTHON} workflow/scripts/40_plot_phytoremediation_framework.py"

rule table1:
    input:
        "data/processed/functional_cohort_638.tsv",
        "results/functional/module28/summary.json",
        "results/functional/module29/species_summary.tsv"
    output:
        "results/tables/module41/table1_genomic_dataset.tsv"
    shell:
        "{PYTHON} workflow/scripts/41_generate_table1.py"

rule table2:
    input:
        "results/functional/module32/gene_variance_partition.tsv"
    output:
        "results/tables/module42/table2_functional_variability.tsv"
    shell:
        "{PYTHON} workflow/scripts/42_generate_table2.py"

rule table3:
    input:
        "results/statistics/module36/marker_associations.tsv"
    output:
        "results/tables/module44/table4_marker_associations.tsv"
    shell:
        "{PYTHON} workflow/scripts/44_generate_table4.py"

rule table4:
    input:
        "results/functional/module28/gene_presence.tsv",
        "results/functional/external_phenotypes/module28/gene_presence.tsv"
    output:
        "results/tables/module43/table3_phenotype_genome_evidence.tsv"
    shell:
        "{PYTHON} workflow/scripts/43_generate_table3.py"


# ============================================================
# Functional analyses using existing genomic matrices
# ============================================================

rule species_distribution:
    input:
        manifest="data/processed/genome_cohort_selected.tsv",
        genes="results/functional/module28/gene_presence.tsv",
        groups="results/functional/module28/group_presence.tsv"
    output:
        prevalence="results/functional/module29/species_gene_prevalence.tsv",
        summary="results/functional/module29/species_summary.tsv"
    shell:
        """
        python workflow/scripts/29_species_functional_distribution.py \
            --manifest {input.manifest} \
            --gene-presence {input.genes} \
            --group-presence {input.groups} \
            --outdir results/functional/module29
        """

rule variance_partition:
    input:
        manifest="data/processed/functional_cohort_638.tsv",
        matrix="results/functional/module28/gene_presence.tsv"
    output:
        "results/functional/module32/gene_variance_partition.tsv"
    shell:
        """
        python workflow/scripts/32_species_variance_partition.py \
            --manifest {input.manifest} \
            --matrix {input.matrix} \
            --outdir results/functional/module32
        """

rule marker_associations:
    input:
        matrix="results/functional/module28/gene_presence.tsv",
        metadata="data/processed/functional_cohort_638.tsv"
    output:
        "results/statistics/module36/marker_associations.tsv"
    params:
        permutations=9999,
        seed=20261026
    shell:
        """
        python workflow/scripts/36_marker_associations.py \
            --matrix {input.matrix} \
            --metadata {input.metadata} \
            --output results/statistics/module36 \
            --permutations {params.permutations} \
            --seed {params.seed}
        """


# ============================================================
# Functional matrix construction (module 28)
# ============================================================

rule functional_matrix:
    input:
        cohort="data/processed/functional_cohort_638_fastas.txt",
        catalog="data/reference/functional_panel/functional_reference_panel_merged.tsv",
        references="data/reference/functional_panel/functional_reference_panel_merged.faa",
        alignments="results/functional/module27/summary.json"
    output:
        genes="results/functional/module28/gene_presence.tsv",
        groups="results/functional/module28/group_presence.tsv",
        summary="results/functional/module28/summary.json"
    params:
        alignments="results/functional/module27",
        max_evalue=1e-5,
        min_qcov=0.7,
        min_pident=0
    shell:
        """
        python workflow/scripts/28_build_functional_matrix.py \
            --alignments {params.alignments} \
            --catalog {input.catalog} \
            --reference-fasta {input.references} \
            --cohort {input.cohort} \
            --outdir results/functional/module28 \
            --max-evalue {params.max_evalue} \
            --min-qcov {params.min_qcov} \
            --min-pident {params.min_pident}
        """


# ============================================================
# Functional homology search (module 27)
# ============================================================

rule functional_homolog_search:
    input:
        cohort="data/processed/functional_cohort_638_fastas.txt",
        references="data/reference/functional_panel/functional_reference_panel_merged.faa"
    output:
        summary="results/functional/module27/summary.json"
    threads: 2
    params:
        sensitivity=5.7,
        evalue="1e-5"
    shell:
        """
        python workflow/scripts/27_search_functional_homologs.py \
            --fasta-list {input.cohort} \
            --references {input.references} \
            --outdir results/functional/module27 \
            --threads {threads} \
            --sensitivity {params.sensitivity} \
            --evalue {params.evalue}
        """


# ============================================================
# Functional cohort validation (module 25)
# ============================================================

rule validate_functional_cohort:
    input:
        cohort="data/processed/genome_cohort_selected.tsv",
        audit="results/ani/module24/representative_proposal_audit.tsv",
        retained="results/ani/module24/proposed_retained.tsv",
        omissions="results/ani/module24/proposed_omissions.tsv",
        pairs="results/ani/module23/ani_pair_observations.tsv",
        phenotype_links="data/processed/phenotype_linked_subset.tsv"
    output:
        manifest="data/processed/functional_cohort_638.tsv",
        fasta_list="data/processed/functional_cohort_638_fastas.txt"
    params:
        genomes_dir="data/raw/genomes/cohort21",
        outdir="results/validation/module25",
        min_ani=99.9,
        min_af_both=90.0,
        expected_original=877,
        expected_retained=638
    shell:
        """
        python workflow/scripts/25_validate_functional_cohort.py \
            --cohort {input.cohort} \
            --audit {input.audit} \
            --retained {input.retained} \
            --omissions {input.omissions} \
            --pairs {input.pairs} \
            --phenotype-links {input.phenotype_links} \
            --genomes-dir {params.genomes_dir} \
            --outdir {params.outdir} \
            --min-ani {params.min_ani} \
            --min-af-both {params.min_af_both} \
            --expected-original {params.expected_original} \
            --expected-retained {params.expected_retained}
        """


rule functional_profiles:
    input:
        groups="results/functional/module28/group_presence.tsv",
        richness="results/functional/module29/genome_functional_richness.tsv"
    output:
        summary="results/functional/module30/summary.json"
    shell:
        """
        python workflow/scripts/30_functional_profiles.py \
            --group-presence {input.groups} \
            --genome-richness {input.richness} \
            --outdir results/functional/module30
        """


rule gene_variability:
    input:
        manifest="data/processed/functional_cohort_638.tsv",
        genes="results/functional/module28/gene_presence.tsv"
    output:
        summary="results/functional/module31/summary.json"
    shell:
        """
        python workflow/scripts/31_gene_variability.py \
            --manifest {input.manifest} \
            --gene-presence {input.genes} \
            --outdir results/functional/module31
        """


rule within_system_codetection:
    input:
        manifest="data/processed/functional_cohort_638.tsv",
        genes="results/functional/module28/gene_presence.tsv"
    output:
        summary="results/functional/module33/summary.json"
    shell:
        """
        python workflow/scripts/33_within_system_codetection.py \
            --gene-presence {input.genes} \
            --manifest {input.manifest} \
            --outdir results/functional/module33
        """


# ============================================================
# Within-species ANI (module 22)
# ============================================================

rule within_species_ani:
    input:
        cohort="data/processed/genome_cohort_selected.tsv"
    output:
        manifest="results/ani/module22/species_ani_manifest.tsv"
    threads: 2
    params:
        fastas="data/raw/genomes/cohort21",
        outdir="results/ani/module22"
    shell:
        """
        python workflow/scripts/22_compute_within_species_ani.py \
            --cohort {input.cohort} \
            --fastas {params.fastas} \
            --output {params.outdir} \
            --threads {threads} \
            --max-species 0
        """

# ============================================================
# ANI redundancy candidates (module 23)
# ============================================================

rule ani_redundancy_candidates:
    input:
        cohort="data/processed/genome_cohort_selected.tsv",
        manifest="results/ani/module22/species_ani_manifest.tsv"
    output:
        groups="results/ani/module23/ani_candidate_members.tsv",
        pairs="results/ani/module23/ani_pair_observations.tsv",
        summary="results/ani/module23/summary.json"
    params:
        sparse_dir="results/ani/module22/sparse",
        outdir="results/ani/module23"
    shell:
        """
        python workflow/scripts/23_ani_redundancy_candidates.py \
            --cohort {input.cohort} \
            --manifest {input.manifest} \
            --sparse-dir {params.sparse_dir} \
            --outdir {params.outdir} \
            --min-ani 99.9 \
            --min-af-both 90
        """

# ============================================================
# ANI representative proposal (module 24)
# ============================================================

rule propose_ani_representatives:
    input:
        cohort="data/processed/genome_cohort_selected.tsv",
        groups="results/ani/module23/ani_candidate_members.tsv",
        pairs="results/ani/module23/ani_pair_observations.tsv",
        phenotype_links="data/processed/phenotype_linked_subset.tsv"
    output:
        retained="results/ani/module24/proposed_retained.tsv",
        omissions="results/ani/module24/proposed_omissions.tsv",
        audit="results/ani/module24/representative_proposal_audit.tsv",
        summary="results/ani/module24/summary.json"
    params:
        outdir="results/ani/module24"
    shell:
        """
        python workflow/scripts/24_propose_ani_representatives.py \
            --cohort {input.cohort} \
            --groups {input.groups} \
            --pairs {input.pairs} \
            --phenotype-links {input.phenotype_links} \
            --outdir {params.outdir} \
            --min-ani 99.9 \
            --min-af-both 90
        """
