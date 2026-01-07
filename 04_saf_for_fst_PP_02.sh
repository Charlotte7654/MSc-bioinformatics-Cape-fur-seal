#!/bin/bash
#PBS -N 04_saf_for_fst_PP_02
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=8:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

# Load ANGSD module
module load app/angsd/0.940

mkdir -p ~/ch2/03_angsd/sfs/

angsd \
  -bam ~/ch2/00_info/bamlist_PP \
  -anc ~/msc/05_refmap/genome_updated/GCA_040869175.1_updated.fasta \
  -out ~/ch2/03_angsd/sfs/03_angsd_saf_PP_02 \
  -doSaf 1 \
  -GL 1 \
  -nThreads 8 \
  -minQ 20 \
  -minMapQ 20 \
  -doCounts 1 \
  -doMajorMinor 1 \
  -doMaf 1
