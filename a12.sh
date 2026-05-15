#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 2
#SBATCH --mem 8G
#SBATCH -J LBAI
#SBATCH --time=0-00:59:59


module load miniconda/24.1.2

source activate /home/user/woon/ML


python update_itp_trustable_dihe_json.py \
  --data ./data \
  --base-itp BK7T.itp \
  --dihe-json dihe_optfile.json \
  --charges type_charge.txt \
  --ff-string ff_string.txt \
  --output-itp BK7T_full_hessfit_scan_refined_relaxed.itp \
  --allow-incomplete \
  --expected-points 8 \
  --allow-high-rmse