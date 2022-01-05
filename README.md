# MPXV_NanoPoreSeq

This is snakefile code accompanying "Nanopore sequencing of a monkeypox virus strain isolated from a pustular lesion in the Central African Republic".


## Preparation of data

```
cd /scratch/
ln ~/Barcode02/barcode02.fastq .
ln ~/NC_003310.fasta .
```

## Execution

```
singularity exec -B $PWD:/data/ -B /bioinfo/:/bioinfo/ -B /scratch/:/scratch/ -B MPXV_NanoPoreSeq:/pipeline/ minion_homopolish_snippy_2_08b.sif  bash -l

snakemake -j10 -s /pipeline/Snakefile3
```
