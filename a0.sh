#!/bin/bash -l
#SBATCH --partition=gpu-k40c
#SBATCH --gres=gpu:2
#SBATCH --job-name=test2
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --time=0-23:59:00
#SBATCH --qos=long

# Define TeraChem variables
export TeraChem=/home/user/woon/terachem-1.95p
export PATH=$TeraChem/bin:$PATH
export LD_LIBRARY_PATH=$TeraChem/lib:$LD_LIBRARY_PATH
export NBOEXE=$TeraChem/bin/nbo6.i4.exe

# Run TeraChem
terachem /home/user/woon/hessfit_work/LBAI/111.ts >111.out

