#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc9
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab21.tpl -e 1PopStab21.est -d -M -n100000 -L20 -c8 -q --seed 243331 -z0.01 -y 3

./fsc28 -t 1PopStab22.tpl -e 1PopStab22.est -d -M -n100000 -L20 -c8 -q --seed 254442 -z0.01 -y 3

./fsc28 -t 1PopStab23.tpl -e 1PopStab23.est -d -M -n100000 -L20 -c8 -q --seed 265553 -z0.01 -y 3

./fsc28 -t 1PopStab24.tpl -e 1PopStab24.est -d -M -n100000 -L20 -c8 -q --seed 276664 -z0.01 -y 3

./fsc28 -t 1PopStab25.tpl -e 1PopStab25.est -d -M -n100000 -L20 -c8 -q --seed 287775 -z0.01 -y 3

