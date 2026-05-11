#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 2
#SBATCH --mem 8G
#SBATCH -J LBAI
#SBATCH --time=0-00:59:59


module load miniconda/24.1.2
module load gaussian/g09

source activate /home/user/woon/ML

export HESSDIR="/home/user/woon/ML/lib/python3.9/site-packages/hessfit"
export HESSPY="/home/user/woon/ML/bin/python"
export PYTHONPATH="$HESSDIR:$PYTHONPATH"
export PATH="$HESSDIR:$PATH"

cd /home/user/woon/hessfit_work/LBAI

"$HESSPY" "$HESSDIR/build_4_hessfit.py" optfile.json --version g09 --path /app/gaussian/g09 --at scratch
