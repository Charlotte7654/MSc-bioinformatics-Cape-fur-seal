#!/bin/bash
#PBS -N 05_LEA_HG_full
#PBS -l select=1:ncpus=16:mem=16GB
#PBS -l walltime=168:00:00
#PBS -m be
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
module load app/R/4.3.2

#in same 00_scripts directory as the script is being submitted 
#I believe cd $PBS_O_WORKDIR means make directory of submission working dir 
R --file=05_LEA_HG_full_Rscript.R
