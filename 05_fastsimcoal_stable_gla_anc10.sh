#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc10
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab26.tpl -e 1PopStab26.est -d -M -n100000 -L20 -c8 -q --seed 298886 -z0.01 -y 3

./fsc28 -t 1PopStab27.tpl -e 1PopStab27.est -d -M -n100000 -L20 -c8 -q --seed 309997 -z0.01 -y 3

./fsc28 -t 1PopStab28.tpl -e 1PopStab28.est -d -M -n100000 -L20 -c8 -q --seed 321108 -z0.01 -y 3

./fsc28 -t 1PopStab29.tpl -e 1PopStab29.est -d -M -n100000 -L20 -c8 -q --seed 332219 -z0.01 -y 3

./fsc28 -t 1PopStab30.tpl -e 1PopStab30.est -d -M -n100000 -L20 -c8 -q --seed 343330 -z0.01 -y 3

