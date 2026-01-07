#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc7
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab11.tpl -e 1PopStab11.est -d -M -n100000 -L20 -c8 -q --seed 132221 -z0.01 -y 3

./fsc28 -t 1PopStab12.tpl -e 1PopStab12.est -d -M -n100000 -L20 -c8 -q --seed 143332 -z0.01 -y 3

./fsc28 -t 1PopStab13.tpl -e 1PopStab13.est -d -M -n100000 -L20 -c8 -q --seed 154443 -z0.01 -y 3

./fsc28 -t 1PopStab14.tpl -e 1PopStab14.est -d -M -n100000 -L20 -c8 -q --seed 165554 -z0.01 -y 3

./fsc28 -t 1PopStab15.tpl -e 1PopStab15.est -d -M -n100000 -L20 -c8 -q --seed 176665 -z0.01 -y 3

