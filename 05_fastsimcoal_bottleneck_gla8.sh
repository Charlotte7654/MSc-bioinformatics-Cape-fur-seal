#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck_gla8
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck_gla

./fsc28 -t 1PopBot16.tpl -e 1PopBot16.est -d -M -n100000 -L20 -c8 -q --seed 187776 -z0.01 -y 3

./fsc28 -t 1PopBot17.tpl -e 1PopBot17.est -d -M -n100000 -L20 -c8 -q --seed 198887 -z0.01 -y 3

./fsc28 -t 1PopBot18.tpl -e 1PopBot18.est -d -M -n100000 -L20 -c8 -q --seed 209998 -z0.01 -y 3

./fsc28 -t 1PopBot19.tpl -e 1PopBot19.est -d -M -n100000 -L20 -c8 -q --seed 221109 -z0.01 -y 3

./fsc28 -t 1PopBot20.tpl -e 1PopBot20.est -d -M -n100000 -L20 -c8 -q --seed 232220 -z0.01 -y 3

