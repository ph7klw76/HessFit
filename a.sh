#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH -J BK7T
#SBATCH --time=2-23:59:59
#SBATCH --qos=long

module load gaussian/g09
source $g09profile


g09 <BK7T.gjf> BK7T.log
formchk -3 BK7T.chk BK7T.fchk



