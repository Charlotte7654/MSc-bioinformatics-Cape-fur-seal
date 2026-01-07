#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck_gla14
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck_gla

./fsc28 -t 1PopBot46.tpl -e 1PopBot46.est -d -M -n100000 -L20 -c8 -q --seed 521106 -z0.01 -y 3

./fsc28 -t 1PopBot47.tpl -e 1PopBot47.est -d -M -n100000 -L20 -c8 -q --seed 532217 -z0.01 -y 3

./fsc28 -t 1PopBot48.tpl -e 1PopBot48.est -d -M -n100000 -L20 -c8 -q --seed 543328 -z0.01 -y 3

./fsc28 -t 1PopBot49.tpl -e 1PopBot49.est -d -M -n100000 -L20 -c8 -q --seed 554439 -z0.01 -y 3

./fsc28 -t 1PopBot50.tpl -e 1PopBot50.est -d -M -n100000 -L20 -c8 -q --seed 565550 -z0.01 -y 3

