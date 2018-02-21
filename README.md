# DMRfinder
A pipeline to identify differentially-methylated regions from whole-genome bisulfite sequencing data

Required arguments:\n
        -n, --chrNum = chromosome number
        -d, --setwd = working directory
        -c, --numCtrl = number of control samples
        -e, --numExp = number of experimental samples
        -a, --genome = genome assembly (hg38, hg19, mm10, rn6, rheMac8)
        -o, --outprefix = title used in output files

Optional arguments:
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

Output arguments:
        --gold_bed = Output bed file of locations for gold DMRs [default = TRUE]
        --silver_bed = Output bed file of locations for silver DMRs [default = TRUE]
        --silver_info = Output txt file of info for silver DMRs [default = TRUE]
        --gold_plots = Output pdf file of plots for gold DMRs [default = TRUE]
        --silver_plots = Output pdf file of plots for silver DMRs [default = FALSE]
        --background = Output bed file of background DMRs [default = FALSE]
        --meth = Output table of smoothed methylation for each sample at each silver DMR [default = FALSE]
        --cov = Output table of total coverage for each sample at each DMR [default = FALSE]
        --background_cov = Output table of total coverage for each sample at each background DMR [default = FALSE]
        --CpGs = Output bed file of tested CpGs [default = FALSE]

Add DSS_file prefixes after all other arguments, with all control samples first (DSS_files/CTRL01_)
