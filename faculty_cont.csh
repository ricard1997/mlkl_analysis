#!/bin/bash
#SBATCH --job-name=pla
#SBATCH -o pla.o%j
#SBATCH -N 1
#SBATCH --ntasks-per-node=56
#SBATCH --cluster=faculty
#SBATCH --partition=vmonje
#SBATCH --qos=vmonje
# SBATCH --nodelist=cpn-v05-[25]
#SBATCH -t 72:00:00            # run time (hh:mm:ss)
#SBATCH --mail-user=ricardox@buffalo.edu
#SBATCH --mail-type=END




source /projects/academic/vmonje/ricardox/.py39/bin/activate

python plantilla.py




