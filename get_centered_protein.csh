#!/bin/bash
#SBATCH --job-name=ggg
#SBATCH -o ggg.o%j
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
# SBATCH --reservation=ubhpc-future
# SBATCH --nodelist=cpn-v05-[25]
#SBATCH -t 24:00:00            # run time (hh:mm:ss)
#SBATCH --mail-user=ricardox@buffalo.edu
#SBATCH --mail-type=END
#SBATCH --constraint=q06|q07|q08|q09

#module load ccrsoft/2023.01
#module load gcc/11.2.0 openmpi/4.1.1
#module load gromacs/2021.5


# SBATCH --job-name=nextraction
# SBATCH -o nextraction.o%j
# SBATCH -N 1
# SBATCH --ntasks-per-node=56
# SBATCH --cluster=faculty
# SBATCH --partition=vmonje
# SBATCH --qos=vmonje
# SBATCH --nodelist=cpn-v05-[25]
# SBATCH -t 72:00:00            # run time (hh:mm:ss)
# SBATCH --mail-user=ricardox@buffalo.edu
# SBATCH --mail-type=END



source /projects/academic/vmonje/ricardox/.py39/bin/activate

python plantilla.py


#echo "SOLU" "SOLU" | gmx trjconv -f prod.xtc -center -s prod.tpr -pbc mol -n index.ndx -o centered_prot.xtc 
#echo "SOLU" | gmx trjconv -f centered_prot.xtc -s prod.tpr -b 0 -e 0 -pbc mol -n index.ndx -o centered_prot.gro
#gmx editconf -f centered_prot.gro -o centered_prot.gro -resnr 1

#echo "Analysis"
#echo "Analysis"
#echo "Analysis"
