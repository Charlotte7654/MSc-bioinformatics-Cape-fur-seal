#!/bin/bash
#PBS -N 03_angsd_SFS_91ind
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=300:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

module load app/angsd/0.940

mkdir -p ~/ch2/03_angsd/sfs

angsd \
  -bam ~/ch2/00_info/bamlist_91_renamed \
  -out ~/ch2/03_angsd/sfs/03_angsd_SFS_91ind \
  -anc ~/msc/05_refmap/genome_updated/GCA_040869175.1_updated.fasta \
  -nThreads 8 \
  -GL 1 \
  -doSaf 1 \
  -minMapQ 20 \
  -minQ 20 \
  -doCounts 1 \
  -doMajorMinor 1 \
  -doMaf 1

realSFS \
-fold 1 \
~/ch2/03_angsd/sfs/03_angsd_SFS_91ind.saf.idx \
-P 8 \
>~/ch2/03_angsd/sfs/03_angsd_SFS_91ind.saf.idx.ml
