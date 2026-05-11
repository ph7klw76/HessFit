#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 50G
#SBATCH -J LBAI-1
#SBATCH --time=2-23:59:59
#SBATCH --qos=long

module load gaussian/g09
source $g09profile


g09 <6_qm.gjf> 6_qm.log
formchk -3 6_qm.chk 6_qm.fchk
g09 <4_qm.gjf> 4_qm.log
formchk -3 4_qm.chk 4_qm.fchk
g09 <5_qm.gjf> 5_qm.log
formchk -3 5_qm.chk 5_qm.fchk





