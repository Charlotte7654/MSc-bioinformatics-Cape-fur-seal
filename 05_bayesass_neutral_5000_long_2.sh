#!/bin/bash

#PBS -N 05_bayesass_neutral_5000_long_2
#PBS -l select=1:ncpus=2:mem=4GB
#PBS -l walltime=744:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
module load app/BayesAss3-SNPs/current

#specified in job submission instead 
#RANDOM_SEED=1234 
TYPE=neutral

#make directory and copy file into it 
mkdir -p ~/ch2/05_bayesass/${TYPE}_5000_long_2/randseed${RANDOM_SEED}/
cp ~/modules/phy2str.pl/ch2.str.immanc ~/ch2/05_bayesass/${TYPE}_5000_long_2/randseed${RANDOM_SEED}/ch2.str.immanc

# Define input and output
INFILE=~/ch2/05_bayesass/${TYPE}_5000_long_2/randseed${RANDOM_SEED}/ch2.str.immanc
RUN_DIR=~/ch2/05_bayesass/${TYPE}_5000_long_2/randseed${RANDOM_SEED}

BA3-SNPS \
--file "$INFILE" \
--loci 5000 \
-o "$RUN_DIR"/randseed${RANDOM_SEED} \
-g \
-b 5000000 \
-i 50000000 \
-n 100 \
-t -u -v \
-m0.25 -a0.60 -f0.025 \
-s "${RANDOM_SEED}"

#unique name for trace and indiv files 
mv "$RUN_DIR"/ch2.str.trace.txt "$RUN_DIR"/randseed${RANDOM_SEED}_5000_long_2.trace.txt
mv "$RUN_DIR"/ch2.str.indiv.txt "$RUN_DIR"/randseed${RANDOM_SEED}_5000_long_2.indiv.txt

