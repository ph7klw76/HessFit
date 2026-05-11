#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 50G
#SBATCH -J LBAI
#SBATCH --time=0-00:59:59

module load gaussian/g09
source $g09profile

grm -f 3_mm_*.log 3_mm_*.chk

for f in 3_mm_*.gjf; do
    base="${f%.gjf}"
    echo "Running $f"
    g09 "$f" "${base}.log"
done