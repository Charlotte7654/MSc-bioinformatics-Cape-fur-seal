#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck_gla12
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck_gla

./fsc28 -t 1PopBot36.tpl -e 1PopBot36.est -d -M -n100000 -L20 -c8 -q --seed 409996 -z0.01 -y 3

./fsc28 -t 1PopBot37.tpl -e 1PopBot37.est -d -M -n100000 -L20 -c8 -q --seed 421107 -z0.01 -y 3

./fsc28 -t 1PopBot38.tpl -e 1PopBot38.est -d -M -n100000 -L20 -c8 -q --seed 432218 -z0.01 -y 3

./fsc28 -t 1PopBot39.tpl -e 1PopBot39.est -d -M -n100000 -L20 -c8 -q --seed 443329 -z0.01 -y 3

./fsc28 -t 1PopBot40.tpl -e 1PopBot40.est -d -M -n100000 -L20 -c8 -q --seed 454440 -z0.01 -y 3

