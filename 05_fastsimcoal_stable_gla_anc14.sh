#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc14
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab46.tpl -e 1PopStab46.est -d -M -n100000 -L20 -c8 -q --seed 521106 -z0.01 -y 3

./fsc28 -t 1PopStab47.tpl -e 1PopStab47.est -d -M -n100000 -L20 -c8 -q --seed 532217 -z0.01 -y 3

./fsc28 -t 1PopStab48.tpl -e 1PopStab48.est -d -M -n100000 -L20 -c8 -q --seed 543328 -z0.01 -y 3

./fsc28 -t 1PopStab49.tpl -e 1PopStab49.est -d -M -n100000 -L20 -c8 -q --seed 554439 -z0.01 -y 3

./fsc28 -t 1PopStab50.tpl -e 1PopStab50.est -d -M -n100000 -L20 -c8 -q --seed 565550 -z0.01 -y 3

