#!/bin/bash
#PBS -N 05_fastsimcoal_bottleneck_gla9
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/bottleneck_gla

./fsc28 -t 1PopBot21.tpl -e 1PopBot21.est -d -M -n100000 -L20 -c8 -q --seed 243331 -z0.01 -y 3

./fsc28 -t 1PopBot22.tpl -e 1PopBot22.est -d -M -n100000 -L20 -c8 -q --seed 254442 -z0.01 -y 3

./fsc28 -t 1PopBot23.tpl -e 1PopBot23.est -d -M -n100000 -L20 -c8 -q --seed 265553 -z0.01 -y 3

./fsc28 -t 1PopBot24.tpl -e 1PopBot24.est -d -M -n100000 -L20 -c8 -q --seed 276664 -z0.01 -y 3

./fsc28 -t 1PopBot25.tpl -e 1PopBot25.est -d -M -n100000 -L20 -c8 -q --seed 287775 -z0.01 -y 3

