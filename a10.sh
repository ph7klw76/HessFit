#!/bin/bash -l
 
#SBATCH -p cpu-epyc-genoa
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 50G
#SBATCH -J LBAI
#SBATCH --time=0-00:59:59

module load gaussian/g09
source $g09profile

g09 < 0_mm_0.gjf > 0_mm_0.log
g09 < 0_mm_192.gjf > 0_mm_192.log
g09 < 0_mm_288.gjf > 0_mm_288.log
g09 < 0_mm_384.gjf > 0_mm_384.log
g09 < 0_mm_480.gjf > 0_mm_480.log
g09 < 0_mm_576.gjf > 0_mm_576.log
g09 < 0_mm_672.gjf > 0_mm_672.log
g09 < 0_mm_768.gjf > 0_mm_768.log
g09 < 0_mm_864.gjf > 0_mm_864.log
g09 < 0_mm_960.gjf > 0_mm_960.log
g09 < 0_mm_96.gjf > 0_mm_96.log
g09 < 1_mm_0.gjf > 1_mm_0.log
g09 < 1_mm_192.gjf > 1_mm_192.log
g09 < 1_mm_288.gjf > 1_mm_288.log
g09 < 1_mm_384.gjf > 1_mm_384.log
g09 < 1_mm_480.gjf > 1_mm_480.log
g09 < 1_mm_576.gjf > 1_mm_576.log
g09 < 1_mm_672.gjf > 1_mm_672.log
g09 < 1_mm_768.gjf > 1_mm_768.log
g09 < 1_mm_864.gjf > 1_mm_864.log
g09 < 1_mm_960.gjf > 1_mm_960.log
g09 < 1_mm_96.gjf > 1_mm_96.log
g09 < 2_mm_0.gjf > 2_mm_0.log
g09 < 2_mm_192.gjf > 2_mm_192.log
g09 < 2_mm_288.gjf > 2_mm_288.log
g09 < 2_mm_384.gjf > 2_mm_384.log
g09 < 2_mm_480.gjf > 2_mm_480.log
g09 < 2_mm_576.gjf > 2_mm_576.log
g09 < 2_mm_672.gjf > 2_mm_672.log
g09 < 2_mm_768.gjf > 2_mm_768.log
g09 < 2_mm_864.gjf > 2_mm_864.log
g09 < 2_mm_960.gjf > 2_mm_960.log
g09 < 2_mm_96.gjf > 2_mm_96.log
g09 < 3_mm_0.gjf > 3_mm_0.log
g09 < 3_mm_192.gjf > 3_mm_192.log
g09 < 3_mm_288.gjf > 3_mm_288.log
g09 < 3_mm_384.gjf > 3_mm_384.log
g09 < 3_mm_480.gjf > 3_mm_480.log
g09 < 3_mm_576.gjf > 3_mm_576.log
g09 < 3_mm_672.gjf > 3_mm_672.log
g09 < 3_mm_768.gjf > 3_mm_768.log
g09 < 3_mm_864.gjf > 3_mm_864.log
g09 < 3_mm_960.gjf > 3_mm_960.log
g09 < 3_mm_96.gjf > 3_mm_96.log
g09 < 4_mm_0.gjf > 4_mm_0.log
g09 < 4_mm_192.gjf > 4_mm_192.log
g09 < 4_mm_288.gjf > 4_mm_288.log
g09 < 4_mm_384.gjf > 4_mm_384.log
g09 < 4_mm_480.gjf > 4_mm_480.log
g09 < 4_mm_576.gjf > 4_mm_576.log
g09 < 4_mm_672.gjf > 4_mm_672.log
g09 < 4_mm_768.gjf > 4_mm_768.log
g09 < 4_mm_864.gjf > 4_mm_864.log
g09 < 4_mm_960.gjf > 4_mm_960.log
g09 < 4_mm_96.gjf > 4_mm_96.log
g09 < 5_mm_0.gjf > 5_mm_0.log
g09 < 5_mm_192.gjf > 5_mm_192.log
g09 < 5_mm_288.gjf > 5_mm_288.log
g09 < 5_mm_384.gjf > 5_mm_384.log
g09 < 5_mm_480.gjf > 5_mm_480.log
g09 < 5_mm_576.gjf > 5_mm_576.log
g09 < 5_mm_672.gjf > 5_mm_672.log
g09 < 5_mm_768.gjf > 5_mm_768.log
g09 < 5_mm_864.gjf > 5_mm_864.log
g09 < 5_mm_960.gjf > 5_mm_960.log
g09 < 5_mm_96.gjf > 5_mm_96.log
g09 < 6_mm_0.gjf > 6_mm_0.log
g09 < 6_mm_192.gjf > 6_mm_192.log
g09 < 6_mm_288.gjf > 6_mm_288.log
g09 < 6_mm_384.gjf > 6_mm_384.log
g09 < 6_mm_480.gjf > 6_mm_480.log
g09 < 6_mm_576.gjf > 6_mm_576.log
g09 < 6_mm_672.gjf > 6_mm_672.log
g09 < 6_mm_768.gjf > 6_mm_768.log
g09 < 6_mm_864.gjf > 6_mm_864.log
g09 < 6_mm_960.gjf > 6_mm_960.log
g09 < 6_mm_96.gjf > 6_mm_96.log
