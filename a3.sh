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


g09 <GauNonBon.gjf> GauNonBon.log
formchk -3 GauNonBon.chk GauNonBon.fchk

g09 <GauHarm.gjf> GauHarm.log
formchk -3 GauHarm.chk GauHarm.fchk




