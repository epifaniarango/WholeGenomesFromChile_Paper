# Average nucleotide diversity within populations
This script was not created by me, but by Roman Briskine
This pipeline calculates average nucleotide diversity for the populations in the South American
expression data set using [pixy](https://pixy.readthedocs.io/en/latest/). The pipeline contains the
following steps.

- Index the input files
- Download data for VQSR
- Call haplotypes for each sample individually
- Genotype all samples jointly
- Perform VQSR
- Apply additional filters (depth, HWE, missingness, MAF)
- Compute nucleotide diversity for sliding windows
- Aggregate the values from all windows


## Input data

- Sequencing data in bam format should be placed into `data/bams`. Index files should also be
  provided. File names should start with the sample name. Sample name is assumed to consist of all
  characters up to the first underscore `_`. Dashes in sample names are allowed.
- Corresponding genome reference `hs37d5.fa.gz` should be placed in `data/ref`
- List of samples should be provided in `data/share/samples.txt`. It should contain a single column
  with sample names. There should be no header line. The pipeline will only process samples specified
  in this file.
- Population mappings should be provided in `data/share/populations.txt`. The file should have no
  headers and provide two columns: sample name and population name.
- Depending on the sequencing depth, different maximum depth thresholds are applied to different
  samples. The threshold is indicated by adding each sample to one of the two files
  `data/share/maxdp53.txt` and `data/share/maxdb62.txt` for 53 and 62 maximum depth thresholds
  respectively.


## How to run


- Create conda environment

```
mamba create -f env/gatk.yml
```

- Optionally create other conda environments used in the pipeline

```
source activate gatk
snakemake --use-conda --conda-create-envs-only
```

- Update the sample list in `data/share/samples.txt` to include all the samples to be processed
- Update the populations file `data/share/populations.txt`
- Create `log` directory with `mkdir -p log`
- Run pipeline. If you are running it on a cluster with SLURM, you can use the provided `run.slurm`
  script. Most likely, you would need to use a slurm
  [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). If your global
  profile is named "slurm", you can schedule the pipeline as

```
sbatch ./run.slurm --profile slurm
```

The output files will be in `data/pi` directory. In particular, the nucleotide diversity estimates
for the fully filtered input data will be in `data/pi/filtered.txt`.

