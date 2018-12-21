#!/bin/bash
#
#SBATCH --workdir=/share/lasallelab/programs/DMRfinder/test
#SBATCH --partition=production			# cluster partition
#SBATCH --mem=48000                   # total memory
#SBATCH --time=1-0							    # time (day-hr)
#SBATCH -n 2								        # cores

module load R

Rscript --vanilla DMRfinder5.R -n 21 -d /share/lasallelab/programs/DMRfinder/test -c 3 -e 3 -a hg38 -o DFtest --mc.cores=6 CTRL1_ CTRL2_ CTRL3_ EXP1_ EXP2_ EXP3_
