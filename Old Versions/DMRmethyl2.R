#!/usr/bin/env Rscript

# Updated 4/20/18
# Authors: Charles Mordaunt
# This code is to be used on cabernet.genomecenter.ucdavis.edu

cat("\n[DMRmethyl]\n\nA pipeline to retrieve smoothed methylation and coverage in predefined differentially-methylated regions 
    from whole-genome bisulfite sequencing data.Note: CpGs are subsetted based on coverage in Ctrl, Exp, and DisCtrl samples.\n\n")

#################################################
# Functions 
#################################################

cat("\n[DMRmethyl] Loading functions\n\n")

print_help_compact <- function (object) {
# Originally written by Trevor Davis as part of optparse package
# Modified by Charles Mordaunt
# Prints help for an Option Parser object, more compact than print_help()
        cat(object@usage, fill = TRUE)
        cat(object@description, fill = TRUE)
        cat("Required arguments:", sep = "\n")
        options_list <- object@options
        for (ii in 1:8) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat("\nOptional arguments:", sep = "\n")
        for (ii in 9:12) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat("\nOutput arguments:", sep = "\n")
        for (ii in 13:14) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat(object@epilogue, fill = TRUE, sep = "")
        return(invisible(NULL))
}

#################################################
# Global Variables
#################################################

# Get arguments from bash script
library(optparse)
cat("[DMRmethyl] Getting arguments from bash script\n\n")
option_list <- list(
        # Required arguments
        make_option(opt_str = c("-n", "--chrNum"), type = "integer", default = NULL, help = "chromosome number"),
        make_option(opt_str = c("-d", "--setwd"), type = "character", default = NULL, help = "working directory"),
        make_option(opt_str = c("-c", "--numCtrl"), type = "integer", default = NULL, help = "number of control samples"),
        make_option(opt_str = c("-e", "--numExp"), type = "integer", default = NULL, help = "number of experimental samples"),
        make_option(opt_str = c("-i", "--numDisCtrl"), type = "integer", default = NULL, help = "number of disease control samples"),
        make_option(opt_str = c("-a", "--genome"), type = "character", default = NULL, help = "genome assembly (hg38, hg19, mm10, rn6, rheMac8)"),
        make_option(opt_str = c("-o", "--outprefix"), type = "character", default = NULL, help = "title used in output files"),
        make_option(opt_str = c("-r", "--DMRfile"), type = "character", default = NULL, help = "bed file of predefined DMRs"),
        
        # Optional arguments
        make_option(opt_str = c("--pctMinCtrl"), type = "double", default = 0.9, help = "minimum percent of control samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--pctMinExp"), type = "double", default = 0.9, help = "minimum percent of experimental samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--pctMinDisCtrl"), type = "double", default = 0.9, help = "minimum percent of disease control samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--mc.cores"), type = "integer", default = 1, help = "cores to use, same as SBATCH -n [default = 1]"),
       
        # Output arguments
        make_option(opt_str = c("--meth"), type = "logical", default = FALSE, help = "Output table of smoothed methylation for each sample at each silver DMR [default = FALSE]"),
        make_option(opt_str = c("--cov"), type = "logical", default = FALSE, help = "Output table of total coverage for each sample at each DMR [default = FALSE]")
)
opt_obj <- OptionParser(option_list = option_list, epilogue = "\nAdd DSS_file prefixes after all other arguments, with all control samples first (DSS_files/CTRL01_)\n", 
                        add_help_option = TRUE, usage = "usage: %prog [arguments]")
opt <- parse_args(object = opt_obj, positional_arguments = c(0, Inf), print_help_and_exit = FALSE)
          
# Test for required arguments
if(opt$options$help){
        print_help_compact(opt_obj)
        stop("", call.=FALSE)
}
if(is.null(opt$options$chrNum)){
        print_help_compact(opt_obj)
        stop("chrNum must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$setwd)){
        print_help_compact(opt_obj)
        stop("setwd must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$numCtrl)){
        print_help_compact(opt_obj)
        stop("numCtrl must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$numExp)){
        print_help_compact(opt_obj)
        stop("numExp must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$genome)){
        print_help_compact(opt_obj)
        stop("genome must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$outprefix)){
        print_help_compact(opt_obj)
        stop("outprefix must be supplied\n\n", call.=FALSE)
}
if(is.null(opt$options$DMRfile)){
        print_help_compact(opt_obj)
        stop("DMRfile must be supplied\n\n", call.=FALSE)
}
if(length(opt$args) < 4){
        print_help_compact(opt_obj)
        stop("At least 4 DSS files must be supplied\n\n", call.=FALSE)
}
cat("\n", str(opt))   

# Assign arguments to global variables
chrNum <- as.numeric(opt$options$chrNum)
setwd(as.character(opt$options$setwd))                    
numCtrl <- as.numeric(opt$options$numCtrl)                  
numExp <- as.numeric(opt$options$numExp) 
numDisCtrl <- as.numeric(opt$options$numDisCtrl)
genome <- as.character(opt$options$genome)
outprefix <- as.character(opt$options$outprefix)
DMRfile <- as.character(opt$options$DMRfile)
numMinCtrl <- ceiling(as.numeric(opt$options$pctMinCtrl)*numCtrl)               
numMinExp <- ceiling(as.numeric(opt$options$pctMinExp)*numExp)
numMinDisCtrl <- ceiling(as.numeric(opt$options$pctMinDisCtrl)*numDisCtrl)
mc.cores <- as.numeric(opt$options$mc.cores)                 
meth <- as.logical(opt$options$meth)
cov <- as.logical(opt$options$cov)
DSSprefix <- opt$args   

################################################
# Packages
################################################
cat("[DMRmethyl] Loading packages\n")
cat("library(BiocGenerics)"); library(BiocGenerics, logical.return = TRUE, quietly = TRUE)
cat("library(Biobase)"); library(Biobase, logical.return = TRUE, quietly = TRUE)
cat("library(S4Vectors)"); library(S4Vectors, logical.return = TRUE, quietly = TRUE)
cat("library(matrixStats)"); library(matrixStats, logical.return = TRUE, quietly = TRUE)
cat("library(DelayedArray)"); library(DelayedArray, logical.return = TRUE, quietly = TRUE)
cat("library(bsseq)"); library(bsseq, logical.return = TRUE, quietly = TRUE)
cat("library(DSS)"); library(DSS, logical.return = TRUE, quietly = TRUE)
cat("library(permute)"); library(permute, logical.return = TRUE, quietly = TRUE)
cat("library(GenomicRanges)"); library(GenomicRanges, logical.return = TRUE, quietly = TRUE)
cat("library(scales)"); library(scales, logical.return = TRUE, quietly = TRUE)

#################################################
# DMRmethyl Pipeline Setup
#################################################

#Set up group variables
CTRLgroup <- paste("C",1:numCtrl,sep="")
EXPgroup <- paste("E",1:numExp,sep="")
DISCTRLgroup <- paste("D",1:numDisCtrl,sep="")

# Load chromosome data
if(genome == "hg38" | genome == "hg19"){chroms = paste("chr",1:22,sep=""); chroms = c(chroms,"chrX","chrY", "chrM")
}else if(genome == "mm10"){chroms = paste("chr",1:19,sep="");chroms = c(chroms,"chrX","chrY", "chrM")  
}else if(genome == "rn6"){chroms = paste("chr",1:20,sep="");chroms = c(chroms,"chrX","chrY", "chrM")
}else if(genome == "rheMac8"){chroms = paste("chr",1:20,sep="");chroms = c(chroms,"chrX","chrY", "chrM")
}else{cat(paste("Warning! Chromosome names unknown because genome is defined as ",genome,"\n"))}

#################################################
# Main DMRmethyl Pipeline
#################################################
# Load DMR file
chrom <- chroms[chrNum]
cat("\n[DMRmethyl] Loading", chrom, "DMR File\n")
DMRs <- read.delim(DMRfile, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(DMRs) <- c("chr", "start", "end")
DMRs$chr <- as.character(DMRs$chr)
DMRs$start <- as.integer(DMRs$start)
DMRs$end <- as.integer(DMRs$end)
DMRs <- subset(DMRs, chr == chrom)

if(dim(DMRs)[1] > 0){
        # Load chromosome DSS files
        cat("\n[DMRmethyl] Loading", chrom, "DSS files\n")
        DSSlist <- list(read.table(paste(DSSprefix[1],chrom,".DSS.txt",sep=""), header=TRUE))
        for(i in 2:length(DSSprefix)) {DSSlist <- c(DSSlist,list(read.table(paste(DSSprefix[i],chrom,".DSS.txt",sep=""), header=TRUE)))}

        # Make BSobject and smooth
        cat("\n[DMRmethyl] Making BSobject and smoothing\n")
        BSobj <- makeBSseqData(DSSlist,c(CTRLgroup,EXPgroup,DISCTRLgroup))
        rm(DSSlist)
        BSobj_smoothed <- BSmooth(BSseq = BSobj, mc.cores = mc.cores, maxGap = 10^8, parallelBy = "sample", ns = 70, h = 1000, mc.preschedule = FALSE, 
                          keep.se = FALSE, verbose = TRUE)    
        cat("\n[DMRmethyl] Completed smoothing\n")

        # Add group assignments to smoothed BSobject
        pData <- pData(BSobj_smoothed)
        pData$type <- as.character(c(rep(c("Ctrl"), numCtrl), rep(c("Exp") ,numExp), rep(c("DisCtrl") ,numDisCtrl)))
        pData(BSobj_smoothed) <- pData

        # Subset BSobject by coverage
        cat("\n[DMRmethyl] Subsetting BSobject by coverage\n")
        BSobj_cov <- getCoverage(BSobj_smoothed)
        keep_loci <- which(rowSums(BSobj_cov[, BSobj_smoothed$type == "Ctrl"] >= 1) >= numMinCtrl 
                           & rowSums(BSobj_cov[, BSobj_smoothed$type == "Exp"] >= 1) >= numMinExp
                           & rowSums(BSobj_cov[, BSobj_smoothed$type == "DisCtrl"] >= 1) >= numMinDisCtrl)
        BSobj_keep <- BSobj_smoothed[keep_loci,]

        # Output smoothed methylation in each DMR for each sample
        if(meth){
                cat("\n[DMRmethyl] Printing predefined DMR smoothed methylation\n")
				if(nrow(DMRs) == 1){
					DMR_meth <- data.frame(t(getMeth(BSseq = BSobj_keep, regions = DMRs, type = "smooth", what = "perRegion"))) # Keep this as a data.frame with 1 row
				} else {
					DMR_meth <- data.frame(getMeth(BSseq = BSobj_keep, regions = DMRs, type = "smooth", what = "perRegion")) 
				}                
				DMR_meth <- round(DMR_meth, 4)
                sample_names <- gsub(".*/", "", DSSprefix)
                sample_names <- gsub("_", "", sample_names)
                colnames(DMR_meth) <- sample_names
                DMR_meth <- cbind(DMRs[,c("chr", "start", "end")], DMR_meth)
                write.table(DMR_meth, paste(outprefix, chrom,"DMR_methylation.txt", sep="_"), sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        }

        # Output coverage in each DMR for each sample
        if(cov){
                cat("\n[DMRmethyl] Printing predefined DMR coverage\n")
				if(nrow(DMRs) == 1){
					DMR_cov <- data.frame(t(getCoverage(BSseq = BSobj_keep, regions = DMRs, type = "Cov", what = "perRegionTotal"))) # Keep this as a data.frame with 1 row
				} else {
					DMR_cov <- data.frame(getCoverage(BSseq = BSobj_keep, regions = DMRs, type = "Cov", what = "perRegionTotal")) 
				}                
				sample_names <- gsub(".*/", "", DSSprefix)
                sample_names <- gsub("_", "", sample_names)
                colnames(DMR_cov) <- sample_names
                DMR_cov <- cbind(DMRs[,c("chr", "start", "end")], DMR_cov)
                write.table(DMR_cov, paste(outprefix, chrom, "DMR_coverage.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
} else{cat("Warning, no DMRs found in ", chrom ,"\n")}

# Cleanup memory
rm(list=ls())
cat("\n[DMRmethyl] Finished!\n\n")
