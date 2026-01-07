#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc12
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab36.tpl -e 1PopStab36.est -d -M -n100000 -L20 -c8 -q --seed 409996 -z0.01 -y 3

./fsc28 -t 1PopStab37.tpl -e 1PopStab37.est -d -M -n100000 -L20 -c8 -q --seed 421107 -z0.01 -y 3

./fsc28 -t 1PopStab38.tpl -e 1PopStab38.est -d -M -n100000 -L20 -c8 -q --seed 432218 -z0.01 -y 3

./fsc28 -t 1PopStab39.tpl -e 1PopStab39.est -d -M -n100000 -L20 -c8 -q --seed 443329 -z0.01 -y 3

./fsc28 -t 1PopStab40.tpl -e 1PopStab40.est -d -M -n100000 -L20 -c8 -q --seed 454440 -z0.01 -y 3

