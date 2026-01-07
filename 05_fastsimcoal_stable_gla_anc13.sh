#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc13
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab41.tpl -e 1PopStab41.est -d -M -n100000 -L20 -c8 -q --seed 465551 -z0.01 -y 3

./fsc28 -t 1PopStab42.tpl -e 1PopStab42.est -d -M -n100000 -L20 -c8 -q --seed 476662 -z0.01 -y 3

./fsc28 -t 1PopStab43.tpl -e 1PopStab43.est -d -M -n100000 -L20 -c8 -q --seed 487773 -z0.01 -y 3

./fsc28 -t 1PopStab44.tpl -e 1PopStab44.est -d -M -n100000 -L20 -c8 -q --seed 498884 -z0.01 -y 3

./fsc28 -t 1PopStab45.tpl -e 1PopStab45.est -d -M -n100000 -L20 -c8 -q --seed 509995 -z0.01 -y 3

