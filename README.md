# MPXV_NanoPoreSeq

This is snakefile code accompanying "Nanopore sequencing of a monkeypox virus strain isolated from a pustular lesion in the Central African Republic".

Vandenbogaert M, Kwasiborski A, Gonofio E, Descorps-Declère S, Selekon B, Nkili Meyong AA, Ouilibona RS, Gessain A, Manuguerra JC, Caro V, Nakoune E, Berthet N. Nanopore sequencing of a monkeypox virus strain isolated from a pustular lesion in the Central African Republic. Sci Rep. 2022 Jun 24;12(1):10768. doi: 10.1038/s41598-022-15073-1. PMID: 35750759; PMCID: PMC9232561.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9232561/

## Installation & Usage

### Installation
Using Docker/Singularity.

All conda/python dependencies are defined in accompanying dependency files:
 - `conda_installed_packages_base.txt`
 - `conda_installed_packages_homopolish.txt`
 - `pip38_installed_packages.txt`

The provided Singularity file is illustrative of the dependency definitions, and on building a target docker/singularity instance.

## Preparation of data

### Basecalling

Input data is supposed to be basecalled, prior to using the provided snakemake file.

Example basecalling instructions (below instructions are uinsg Guppy v 3.2.4, and are indicative only):

Example using CPUs:

```
dir=/opt/Guppy/ont-guppy-cpu_3.4.4/ont-guppy-cpu/bin

${dir}/guppy_basecaller --kit ${kit} --flowcell ${flowcell} --barcode_kits ${barcode_kit} -i ${indir}/ -s ${outdir} --num_callers 4 --cpu_threads_per_caller 20 -q 4000 --qscore_filtering --min_qscore ${min_qscore} --disable_pings --trim_barcodes

```

Example using GPUs:

Works on Tesla P100 only.

```
${dir}/guppy_basecaller -i /data/fast5_pass/ --save_path /scratch/out/ --flowcell ${flowcell} --kit ${barcode_kit} --gpu_runners_per_device 8 -r --qscore_filtering --min_qscore 7 -x auto --disable_pings --trim_barcodes

```

### Organization of FASTQ files
and reference genome (here reference NC_003310).

Working directory will be `/scratch/`.

```
cd /scratch/
ln ~/RawData/*.fastq .
ln ~/NC_003310.fasta .
```

### Execution

```
snakemake -j10 -s Snakefile
```

Jupyter python/R notebooks for downstream analysis will be added shortly.
