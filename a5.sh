#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 50G
#SBATCH -J LBAI
#SBATCH --time=2-23:59:59
#SBATCH --qos=long

module load gaussian/g09
source $g09profile


g09 <hessfit4gau.gjf> hessfit4gau.log
formchk -3 hessfit4gau.chk hessfit4gau.fchk




