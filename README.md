# DMRfinder

#### A pipeline to identify differentially-methylated regions from whole-genome bisulfite sequencing data

By Charles Mordaunt

Readme updated 7/12/18

## Installing required packages:
First, start a new R session and install the packages below.

        install.packages(c("optparse","permute","scales","MASS","mgcv","survival"))
        source("https://bioconductor.org/biocLite.R")
        biocLite()
        biocLite(c("BiocGenerics","Biobase","S4Vectors","matrixStats","DelayedArray","bsseq","DSS",
                   "GenomicRanges"))

## DSS_file format
Input files should be split by each sample and each chromosome. The format should be in DSS which is a tab-delimited text file with 4 columns: chromosome (chr), position (pos), total reads (N), and methylated reads (X) (see below).
        
        chr	        pos	N	X
        chr21	5013971	1	1
        chr21	5014046	1	1
        chr21	5014056	1	1
        chr21	5014082	1	1
        chr21	5014097	1	0


## Required arguments:  
These arguments are required for running DMRfinder

        -n, --chrNum = chromosome number (1-25, depending on genome, last 3 are chrX, chrY, chrM) 
        -d, --setwd = working directory  
        -c, --numCtrl = number of control samples  
        -e, --numExp = number of experimental samples  
        -a, --genome = genome assembly (hg38, hg19, mm10, rn6, rheMac8)  
        -o, --outprefix = title used in output files  

## Optional arguments:  
These arguments generally don't need to be changed. For parallel processing for smoothing, increase mc.cores. For more stringent DMR criteria, increase pctMinCtrl and pctMinExp to 1 (CpGs have 1 read in all samples) and increase meanDiff_cutoff to 0.1 (DMRs have at least 10% difference between groups).

        --pctMinCtrl = minimum percent of control samples with 1 read at CpG [default = 0.9]  
        --pctMinExp = minimum percent of experimental samples with 1 read at CpG [default = 0.9]  
        --mc.cores = cores to use, same as SBATCH -n [default = 1]  
        --estimate.var = method to estimate variance for t-test (same, paired, group2) [default = same]  
        --meanDiff_cutoff = minimum difference between group means for DMRs [default = 0.05]  
        --maxGap = maximum distance between all consecutive CpGs in a DMR [default = 300]  
        --invdensity_cutoff = maximum average distance between consecutive CpGs in a DMR [default = 300]  
        --colorCtrl = color for control samples in plots [default = 3366CC]  
        --colorExp = color for experimental samples in plots [default = FF3366]  
        --nperm = number of permutations to perform for FWER estimation [default = 1000]  
  
## Output arguments:  
These arguments control the output files, including the DMR locations, plots, methylation, and coverage.

        --gold_bed = Output bed file of locations for gold DMRs [default = TRUE]  
        --silver_bed = Output bed file of locations for silver DMRs [default = TRUE]  
        --silver_info = Output txt file of info for silver DMRs [default = TRUE]  
        --gold_plots = Output pdf file of plots for gold DMRs [default = TRUE]  
        --silver_plots = Output pdf file of plots for silver DMRs [default = FALSE]  
        --background = Output bed file of background DMRs [default = FALSE]  
        --meth = Output table of smoothed methylation for each sample at each silver DMR [default = FALSE]  
        --cov = Output table of total coverage for each sample at each DMR [default = FALSE]  
        --background_cov = Output table of total coverage for each sample at each background DMR
          [default = FALSE]  
        --CpGs = Output bed file of tested CpGs [default = FALSE]  
  
#### Add DSS_file prefixes after all other arguments, with all control samples first (DSS_files/CTRL01_)  

## Example call:
This call is submitted from the shell and runs DMRfinder on chromosome 21, comparing 3 control and 3 experimental umbilical cord blood samples.

        Rscript --vanilla /share/lasallelab/programs/DMRfinder/DMRfinder5.R -n 21 -d /share/lasallelab/Charles/CM_WGBS_CordBlood -c 3 -e 3 -a hg38 -o DMRfinder DSS_files/JLCM001A_ DSS_files/JLCM001C_ DSS_files/JLCM002A_ DSS_files/JLCM001B_ DSS_files/JLCM001D_ DSS_files/JLCM002B_

![Example DMR plot](https://github.com/cemordaunt/DMRfinder/blob/master/DMRplot.png)
Example plot of DMR methylation
