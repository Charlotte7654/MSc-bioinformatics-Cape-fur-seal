#!/bin/bash
#PBS -N 05_fastsimcoal_stable_gla_anc6
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=24:00:00
#PBS -m e
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
cd ~/software/fastsimcoal/stable_gla_anc

./fsc28 -t 1PopStab6.tpl -e 1PopStab6.est -d -M -n100000 -L20 -c8 -q --seed 76666 -z0.01 -y 3

./fsc28 -t 1PopStab7.tpl -e 1PopStab7.est -d -M -n100000 -L20 -c8 -q --seed 87777 -z0.01 -y 3

./fsc28 -t 1PopStab8.tpl -e 1PopStab8.est -d -M -n100000 -L20 -c8 -q --seed 98888 -z0.01 -y 3

./fsc28 -t 1PopStab9.tpl -e 1PopStab9.est -d -M -n100000 -L20 -c8 -q --seed 109999 -z0.01 -y 3

./fsc28 -t 1PopStab10.tpl -e 1PopStab10.est -d -M -n100000 -L20 -c8 -q --seed 121110 -z0.01 -y 3

