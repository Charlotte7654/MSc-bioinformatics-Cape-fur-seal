#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck13
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck

./fsc28 -t 1PopBot41.tpl -e 1PopBot41.est -d -M -n100000 -L20 -c8 -q --seed 465551 -z0.01 -y 3

./fsc28 -t 1PopBot42.tpl -e 1PopBot42.est -d -M -n100000 -L20 -c8 -q --seed 476662 -z0.01 -y 3

./fsc28 -t 1PopBot43.tpl -e 1PopBot43.est -d -M -n100000 -L20 -c8 -q --seed 487773 -z0.01 -y 3

./fsc28 -t 1PopBot44.tpl -e 1PopBot44.est -d -M -n100000 -L20 -c8 -q --seed 498884 -z0.01 -y 3

./fsc28 -t 1PopBot45.tpl -e 1PopBot45.est -d -M -n100000 -L20 -c8 -q --seed 509995 -z0.01 -y 3

