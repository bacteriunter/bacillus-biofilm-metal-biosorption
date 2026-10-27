#!/usr/bin/env bash

OUT="software_versions.txt"

{
    echo "=== SYSTEM ==="
    date -Is
    uname -a
    cat /etc/os-release
    echo

    echo "=== ENVIRONMENT ==="
    echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-not_active}"
    echo "CONDA_PREFIX=${CONDA_PREFIX:-not_active}"
    echo "PATH=$PATH"
    echo

    echo "=== SOFTWARE VERSIONS ==="

    check_version() {
        local name="$1"
        shift
        echo "--- $name ---"
        if command -v "$1" >/dev/null 2>&1; then
            "$@" 2>&1 | head -n 8
        else
            echo "NOT FOUND IN CURRENT PATH"
        fi
        echo
    }

    check_version "Python" python --version
    check_version "Python3" python3 --version
    check_version "Conda" conda --version
    check_version "Snakemake" snakemake --version
    check_version "NCBI datasets" datasets version
    check_version "NCBI dataformat" dataformat version
    check_version "SeqKit" seqkit version
    check_version "skani" skani --version
    check_version "MMseqs2" mmseqs version
    check_version "CheckM" checkm version
    check_version "CheckM2" checkm2 --version
    check_version "Prodigal" prodigal -v
    check_version "Prokka" prokka --version
    check_version "Bakta" bakta --version
    check_version "BLAST+" blastp -version
    check_version "HMMER" hmmsearch -h
    check_version "GNU parallel" parallel --version
    check_version "R" R --version
    check_version "Git" git --version

    echo "=== PYTHON PACKAGES ==="
    python3 - <<'PY'
from importlib import metadata

packages = [
    "numpy",
    "pandas",
    "scipy",
    "biopython",
    "scikit-learn",
    "statsmodels",
    "matplotlib",
    "seaborn",
    "snakemake",
    "requests",
    "tqdm",
    "pyyaml",
    "networkx",
    "openpyxl",
]

for name in packages:
    try:
        print(f"{name}: {metadata.version(name)}")
    except metadata.PackageNotFoundError:
        print(f"{name}: NOT INSTALLED IN CURRENT PYTHON")
PY

    echo
    echo "=== R PACKAGES ==="
    if command -v Rscript >/dev/null 2>&1; then
        Rscript -e '
        pkgs <- c("ggplot2","dplyr","tidyr","readr",
                  "vegan","ape","phytools","pheatmap")
        for (p in pkgs) {
          v <- if (requireNamespace(p, quietly=TRUE))
            as.character(packageVersion(p)) else "NOT INSTALLED"
          cat(p, ":", v, "\n")
        }'
    else
        echo "Rscript NOT FOUND"
    fi

    echo
    echo "=== CONDA ENVIRONMENTS ==="
    if command -v conda >/dev/null 2>&1; then
        conda env list
    fi

    echo
    echo "=== PROJECT FILES ==="
    find workflow -maxdepth 3 -type f \
        \( -name '*.py' -o -name '*.sh' -o \
           -name 'Snakefile' -o -name '*.yaml' \
           -o -name '*.yml' \) \
        2>/dev/null | sort

} > "$OUT" 2>&1

echo "Versiones guardadas en: $(pwd)/$OUT"
