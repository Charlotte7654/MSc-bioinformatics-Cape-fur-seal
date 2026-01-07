#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck_gla_anc11
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck_gla_anc

./fsc28 -t 1PopBot31.tpl -e 1PopBot31.est -d -M -n100000 -L20 -c8 -q --seed 354441 -z0.01 -y 3

./fsc28 -t 1PopBot32.tpl -e 1PopBot32.est -d -M -n100000 -L20 -c8 -q --seed 365552 -z0.01 -y 3

./fsc28 -t 1PopBot33.tpl -e 1PopBot33.est -d -M -n100000 -L20 -c8 -q --seed 376663 -z0.01 -y 3

./fsc28 -t 1PopBot34.tpl -e 1PopBot34.est -d -M -n100000 -L20 -c8 -q --seed 387774 -z0.01 -y 3

./fsc28 -t 1PopBot35.tpl -e 1PopBot35.est -d -M -n100000 -L20 -c8 -q --seed 398885 -z0.01 -y 3

