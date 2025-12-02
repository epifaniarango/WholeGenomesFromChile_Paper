#!/usr/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=7-00:00:00
#SBATCH --output=log/main_%j

module load mamba

eval "$(conda shell.bash hook)"
conda activate gatk

snakemake --use-conda "$@"
