# DMRfinder

#### A pipeline to identify differentially-methylated regions from whole-genome bisulfite sequencing data

By Charles Mordaunt

Readme updated 12/20/18

## Running DMRfinder:
First prepare WGBS data in DSS file format for at least two groups of two samples. Then run DMRfinder.R with the appropriate arguments as in DMRfinder_example.sh. As part of the pipeline, DMRfinder.R will check and install needed packages.

## DSS file format:
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
        --nperm = number of permutations to perform for FWER estimation [default = 1000]  
  
## Output arguments:  
These arguments control the output files, including the DMR locations, plots, methylation, and coverage.

        --silver_bed = Output bed file of locations for silver DMRs [default = TRUE]  
        --silver_info = Output txt file of info for silver DMRs [default = TRUE]  
        --gold_plots = Output pdf file of plots for gold DMRs [default = TRUE]
        --background = Output bed file of background DMRs [default = TRUE]  
        --gold_bed = Output bed file of locations for gold DMRs [default = TRUE]  
        --silver_plots = Output pdf file of plots for silver DMRs [default = FALSE]  
        --meth = Output table of smoothed methylation for each sample at each silver DMR 
          [default = FALSE]  
        --cov = Output table of total coverage for each sample at each DMR [default = FALSE]  
        --background_cov = Output table of total coverage for each sample at each background DMR
          [default = FALSE]  
        --CpGs = Output bed file of tested CpGs [default = FALSE]  
  
#### Add DSS_file prefixes after all other arguments, with all control samples first (DSS_files/CTRL01_)  

## Example call:
This call is submitted from the shell and runs DMRfinder on chromosome 21, comparing 3 control and 3 experimental samples.

        Rscript --vanilla DMRfinder.R -n 21 -d /share/lasallelab/programs/DMRfinder/test -c 3 -e 3 
        -a hg38 -o DFtest CTRL1_ CTRL2_ CTRL3_ EXP1_ EXP2_ EXP3_

## Example output: 
### DMR info

chr | start | end | CpGs | width | invdensity | areaStat | maxStat | tstat_sd | Ctrl_mean | Exp_mean | meanDiff | direction | FWER | Rel_FWER
---|---|---|---|---|---|---|---|---|---|---|---|---|---|---
chr21 | 46699259 | 46699321 | 16 | 62 | 3.88 | -143.18 | -8.76 | 0.07 | 0.73 | 0.16 | -0.58 | hypo | 0 | 0.00
chr21 | 43779213 | 43779862 | 9 | 649 | 72.11 | 42.30 | 5.35 | 0.06 | 0.44 | 0.75 | 0.31 | hyper | 14 | 0.78
chr21 | 45482332 | 45482675 | 9 | 343 | 38.11 | 29.90 | 3.50 | 0.06 | 0.52 | 0.78 | 0.26 | hyper | 18 | 1.00
chr21 | 34890783 | 34890828 | 10 | 45 | 4.50 | 28.18 | 2.85 | 0.06 | 0.13 | 0.29 | 0.16 | hyper | 18 | 1.00
chr21 | 43985999 | 43986281 | 8 | 282 | 35.25 | 26.54 | 4.15 | 0.06 | 0.16 | 0.35 | 0.20 | hyper | 18 | 1.00

### DMR plot
![Example DMR plot](https://github.com/cemordaunt/DMRfinder/blob/master/DMRplot.png)

## Session Info:
        R version 3.5.0 (2018-04-23)
        Platform: x86_64-pc-linux-gnu (64-bit)
        Running under: Ubuntu 16.04.5 LTS

        Matrix products: default
        BLAS/LAPACK: /usr/lib/libopenblasp-r0.2.18.so

        locale:
        [1] LC_CTYPE=en_US       LC_NUMERIC=C         LC_TIME=en_US       
        [4] LC_COLLATE=en_US     LC_MONETARY=en_US    LC_MESSAGES=en_US   
        [7] LC_PAPER=en_US       LC_NAME=C            LC_ADDRESS=C        
        [10] LC_TELEPHONE=C       LC_MEASUREMENT=en_US LC_IDENTIFICATION=C 

        attached base packages:
        [1] splines   stats4    parallel  stats     graphics  grDevices utils    
        [8] datasets  methods   base     

        other attached packages:
        [1] scales_1.0.0                permute_0.9-4              
        [3] DSS_2.30.0                  bsseq_1.18.0               
        [5] SummarizedExperiment_1.12.0 GenomicRanges_1.34.0       
        [7] GenomeInfoDb_1.18.1         DelayedArray_0.8.0         
        [9] BiocParallel_1.16.2         IRanges_2.16.0             
        [11] matrixStats_0.54.0          S4Vectors_0.20.1           
        [13] Biobase_2.42.0              BiocGenerics_0.28.0        
        [15] magrittr_1.5                remotes_2.0.2              
        [17] BiocManager_1.30.4         

        loaded via a namespace (and not attached):
        [1] Rcpp_1.0.0               compiler_3.5.0           XVector_0.22.0          
        [4] R.methodsS3_1.7.1        R.utils_2.7.0            bitops_1.0-6            
        [7] tools_3.5.0              DelayedMatrixStats_1.4.0 zlibbioc_1.28.0         
        [10] rhdf5_2.26.1             lattice_0.20-35          BSgenome_1.50.0         
        [13] Matrix_1.2-14            GenomeInfoDbData_1.2.0   rtracklayer_1.42.1      
        [16] Biostrings_2.50.1        gtools_3.8.1             locfit_1.5-9.1          
        [19] grid_3.5.0               data.table_1.11.8        HDF5Array_1.10.1        
        [22] XML_3.98-1.16            limma_3.38.3             Rhdf5lib_1.4.2          
        [25] GenomicAlignments_1.18.0 Rsamtools_1.34.0         colorspace_1.3-2        
        [28] RCurl_1.95-4.11          munsell_0.5.0            R.oo_1.22.0 
