#!/bin/bash
#PBS -N 03.64_angsd_GL
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=4:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

module load app/angsd/0.940

mkdir -p ~/ch2/03_angsd/

angsd \
  -bam ~/ch2/00_info/bamlist_93_renamed \
  -GL 1 \
  -ref ~/msc/05_refmap/genome_updated/GCA_040869175.1_updated.fasta \
  -out ~/ch2/03_angsd/03.64_angsd_GL \
  -nThreads 8 \
  -doCounts 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doPost 1 \
  -doBcf 1 \
  -minQ 20 \
  -minMapQ 20 \
  -SNP_pval 1e-4 \
  -minInd 47 \
  -minMaf 0.05 \
  -doSNPstat 1 \
  -doHWE 1 \
  --ignore-RG 0 \
  -doDepth 1 \
  -dumpCounts 2

