#!/usr/bin/env Rscript

# Updated 5/31/17
# Authors: Keith Dunaway and Charles Mordaunt
# This code is to be used on cabernet.genomecenter.ucdavis.edu

cat("\n[DMRplot]\n\nPlot smoothed methylation in target regions from whole-genome bisulfite sequencing data\n\n")

#################################################
# Functions 
#################################################

# Functions ###
cat("\n[DMRplot] Loading functions\n\n")

print_help_compact <- function (object) {
        # Originally written by Trevor Davis as part of optparse package
        # Modified by Charles Mordaunt
        # Prints help for an Option Parser object, more compact than print_help()
        cat(object@usage, fill = TRUE)
        cat(object@description, fill = TRUE)
        cat("Required arguments:", sep = "\n")
        options_list <- object@options
        for (ii in 1:7) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat("\nOptional arguments:", sep = "\n")
        for (ii in 8:17) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat(object@epilogue, fill = TRUE, sep = "")
        return(invisible(NULL))
}

# Function written by Dave Tang ###
# Downloaded from: https://github.com/davetang/bedr/blob/master/R/bed_to_granges.R

# Converts bed file to GRanges object
bed_to_granges <- function(file){
        df <- read.table(file, header=F, stringsAsFactors=F)
        if(length(df) > 6){
                df <- df[,-c(7:length(df))]
        }
        if(length(df)<3){
                stop("File has less than 3 columns")
        }
        header <- c('chr','start','end','id','score','strand')
        names(df) <- header[1:length(names(df))]
        if('strand' %in% colnames(df)){
                df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
        }
        if(length(df)==3){
                gr <- with(df, GRanges(chr, IRanges(start, end)))
        } else if (length(df)==4){
                gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
        } else if (length(df)==5){
                gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
        } else if (length(df)==6){
                gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
        }
        return(gr)
}

# Functions originally by Kasper Daniel Hansen, modified by Charles Mordaunt
plotGeneTrack2 <- function (gr, geneTrack, cex) 
{
        geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
        ol <- findOverlaps(geneTrack_gr, gr)
        genes <- geneTrack[queryHits(ol), ]
        plot(start(gr), rep(1, length(start(gr))), type = "n", xaxt = "n", yaxt = "n", bty = "n",   
             ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)), xlab = "",          
             ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
        mtext("Genes", side = 2, at = 1.25, las = 1, line = 1, cex = 0.8, adj = 0.75) 
        mtext("+", side = 2, at = 0.65, las = 1, line = 1, cex = 1, adj = 0)
        mtext("-", side = 2, at = -0.45, las = 1, line = 1, cex = 1, adj = 0)
        if (nrow(genes) > 0) {
                for (g in 1:nrow(genes)) {
                        geneind2 = which(geneTrack$gene_name == genes$gene_name[g])
                        geneind2 = min(geneind2)
                        direction = unique(geneTrack$strand[geneind2])
                        ES = geneTrack$start[geneind2]
                        EE = geneTrack$end[geneind2]
                        Exons = cbind(ES, EE)
                        if (direction == "+") {
                                apply(Exons, 1, function(x) polygon(c(x[1], x[2], x[2], x[1]), c(0.45, 0.45, 0.85, 0.85), col = "#000099"))
                                text((x = max(start(gr), min(ES)) + min(end(gr), max(EE)))/2, y = 1.2, labels = genes$gene_name[g], cex = cex, font = 3)
                        }
                        else {
                                apply(Exons, 1, function(x) polygon(c(x[1], x[2], x[2], x[1]), c(-0.25, -0.25, -0.65, -0.65), col = "#000099"))
                                text(x = (max(start(gr), min(ES)) + min(end(gr), max(EE)))/2, y = -1, labels = genes$gene_name[g], cex = cex, font = 3)
                        }
                }
        }
}

plotRegion2 <- function (BSseq, region = NULL, extend = 0, main = "", addRegions = NULL, 
          annoTrack = NULL, cex.anno = 1, geneTrack = NULL, cex.gene = 1.5, 
          col = NULL, lty = NULL, lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected", 
          stat.col = "black", stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8, 8), mainWithWidth = TRUE, regionCol = alpha("red", 0.1), 
          addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE) 
{
        opar <- par(mar = c(0, 4.1, 0, 0), oma = c(2, 1, 2, 2), mfrow = c(1,1))
        on.exit(par(opar))
        if (is.null(BSseqStat) && !is.null(annoTrack) && !is.null(geneTrack)){layout(matrix(1:3, ncol = 1), heights = c(2, 0.5, 2))}
        if (is.null(annoTrack) && !is.null(BSseqStat) && !is.null(geneTrack)){layout(matrix(1:3, ncol = 1), heights = c(2, 2, 2))}
        if (is.null(geneTrack) && !is.null(BSseqStat) && !is.null(annoTrack)){layout(matrix(1:3, ncol = 1), heights = c(2, 2, 0.5))}
        if (is.null(BSseqStat) && is.null(annoTrack) && !is.null(geneTrack)){layout(matrix(1:2, ncol = 1), heights = c(2, 2))}
        if (is.null(BSseqStat) && is.null(geneTrack) && !is.null(annoTrack)){layout(matrix(1:2, ncol = 1), heights = c(2, 0.5))}
        if (is.null(annoTrack) && is.null(geneTrack) && !is.null(BSseqStat)){layout(matrix(1:2, ncol = 1), heights = c(2, 2))}
        if (is.null(annoTrack) && is.null(geneTrack) && is.null(BSseqStat)){layout(matrix(1:1, ncol = 1), heights = c(2))}
        else {layout(matrix(1:4, ncol = 1), heights = c(2, 2, 0.5, 2))}
        .plotSmoothData2(BSseq = BSseq, region = region, extend = extend, 
                        addRegions = addRegions, col = col, lty = lty, lwd = lwd, 
                        regionCol = regionCol, addTicks = addTicks, addPoints = addPoints, 
                        pointsMinCov = pointsMinCov, highlightMain = highlightMain)
        gr <- bsseq:::.bsGetGr(BSseq, region, extend)
        if (!is.null(BSseqStat)) {
                if (is(BSseqStat, "BSseqTstat")) {
                        stat.values <- as.array(getStats(BSseqStat)[, "tstat.corrected"])
                        stat.values <- as.array(stat.values)
                        stat.type <- "tstat"
                }
                if (is(BSseqStat, "BSseqStat")) {
                        stat.type <- getStats(BSseqStat, what = "stat.type")
                        if (stat.type == "tstat") {
                                stat.values <- getStats(BSseqStat, what = "stat")
                                stat.values <- as.array(stat.values)
                        }
                        if (stat.type == "fstat") {
                                stat.values <- sqrt(getStats(BSseqStat, what = "stat"))
                                stat.values <- as.array(stat.values)
                        }
                }
                plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n", 
                     ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "", 
                     ylab = stat.type, cex.lab = 1.3)
                axis(side = 2, at = c(-5, 0, 5), cex.axis = 1.2)
                abline(h = 0, col = "grey60")
                bsseq:::.bsPlotLines(start(BSseqStat), stat.values, lty = stat.lty, 
                             col = stat.col, lwd = stat.lwd, plotRange = c(start(gr), end(gr)))
                bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(-5.75, 5.75), regionCol = regionCol, highlightMain = highlightMain)
        }
        if (!is.null(annoTrack)) 
                plotAnnoTrack2(gr, annoTrack, cex.anno)
        if (!is.null(geneTrack)) 
                plotGeneTrack2(gr, geneTrack, cex.gene)
        if (!is.null(main)) {
                main <- bsseq:::.bsPlotTitle(gr = region, extend = extend, main = main, mainWithWidth = mainWithWidth)
                mtext(side = 3, text = main, outer = TRUE, cex = 0.9, padj = 0)
        }
        return(invisible(NULL))
}

plotManyRegions2 <- function (BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL, 
          annoTrack = NULL, cex.anno = 1, geneTrack = NULL, cex.gene = 1.5, 
          col = NULL, lty = NULL, lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected", 
          stat.col = "black", stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8, 8), mainWithWidth = TRUE, regionCol = alpha("red", 0.1), 
          addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE, verbose = TRUE) 
{
        cat("[plotManyRegions] preprocessing ...")
        if (!is.null(regions)) {
                if (is(regions, "data.frame")) 
                        gr <- data.frame2GRanges(regions, keepColumns = FALSE)
                else gr <- regions
                if (!is(gr, "GRanges")) 
                        stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
        }
        else {
                gr <- granges(BSseq)
        }
        gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
        BSseq <- subsetByOverlaps(BSseq, gr)
        if (length(start(BSseq)) == 0) 
                stop("No overlap between BSseq data and regions")
        if (!is.null(main) && length(main) != length(gr)) 
                main <- rep(main, length = length(gr))
        cat("done\n")
        for (ii in seq(along = gr)) {
                if (verbose) 
                        cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n", 
                                    ii, nrow(regions)))
                plotRegion2(BSseq = BSseq, region = regions[ii, ], extend = extend, 
                           col = col, lty = lty, lwd = lwd, main = main[ii], 
                           BSseqStat = BSseqStat, stat = stat, stat.col = stat.col, 
                           stat.lwd = stat.lwd, stat.lty = stat.lty, stat.ylim = stat.ylim, 
                           addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth, 
                           annoTrack = annoTrack, cex.anno = cex.anno, geneTrack = geneTrack, 
                           cex.gene = cex.gene, addTicks = addTicks, addPoints = addPoints, 
                           pointsMinCov = pointsMinCov, highlightMain = highlightMain)
        }
}

.plotSmoothData2 <- function (BSseq, region, extend, addRegions, col, lty, lwd, regionCol, 
          addTicks, addPoints, pointsMinCov, highlightMain) 
{
        gr <- bsseq:::.bsGetGr(BSseq, region, extend)
        BSseq <- subsetByOverlaps(BSseq, gr)
        sampleNames <- sampleNames(BSseq)
        names(sampleNames) <- sampleNames
        positions <- start(BSseq)
        smoothPs <- getMeth(BSseq, type = "smooth")
        rawPs <- getMeth(BSseq, type = "raw")
        coverage <- getCoverage(BSseq)
        if (addPoints) {
                rawPs <- as.array(rawPs)
                coverage <- as.array(coverage)
        }
        smoothPs <- as.array(smoothPs)
        ymin <- if(floor(min(smoothPs)*10)/10 < 0){0} else {floor(min(smoothPs)*10)/10}
        ymax <- if(ceiling(max(smoothPs)*10)/10 > 1){1} else {ceiling(max(smoothPs)*10)/10}
        ymid <- round(mean(c(ymin, ymax)), 2)
        colEtc <- bsseq:::.bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)
        plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n", 
             ylim = c(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.05), xlim = c(start(gr), end(gr)), xlab = "", 
             ylab = "Methylation", cex.lab = 1.3)
        axis(side = 2, at = c(ymin, ymid, ymax), cex.axis = 1.2, labels = c(ymin*100, ymid*100, ymax*100))
        if (addTicks) 
                rug(positions, ticksize = 0.08)
        bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0, 1), regionCol = regionCol, highlightMain = highlightMain)
        if (addPoints) {
                sapply(1:ncol(BSseq), function(sampIdx) {
                        abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
                })
        }
        sapply(1:ncol(BSseq), function(sampIdx) {  #Add alpha to col
                bsseq:::.bsPlotLines(positions, smoothPs[, sampIdx], col = alpha(colEtc$col[sampIdx], 0.7), 
                             lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx], 
                             plotRange = c(start(gr), end(gr)))
        })
        if (addPoints) {
                sapply(1:ncol(BSseq), function(sampIdx) {
                        bsseq:::.bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx], col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov)
                })
        }
}

plotAnnoTrack2 <- function (gr, annoTrack, cex) # make thicker, add border
{
        if (!all(sapply(annoTrack, function(xx) is(xx, "GRanges")))) 
                stop("all elements in 'annoTrack' needs to be 'GRanges'")
        plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
             ylim = c(0.5, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), xlab = "", ylab = "")
        lapply(seq(along = annoTrack), function(ii) {
                jj <- length(annoTrack) + 1 - ii
                ir <- subsetByOverlaps(annoTrack[[ii]], gr)
                if (length(ir) > 0) 
                        rect(start(ir) - 0.5, jj - 0.3, end(ir), jj + 0.3, col = "#006600", lwd = 0.8)
                mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1, cex = cex, adj = 0.6)
        })
}
#################################################
# Global Variables
#################################################

# Get arguments from bash script
library(optparse)
cat("[DMRplot] Getting arguments from bash script\n\n")
option_list <- list(
        # Required arguments
        make_option(opt_str = c("-n", "--chrNum"), type = "integer", default = NULL, help = "chromosome number [required]"),
        make_option(opt_str = c("-d", "--setwd"), type = "character", default = NULL, help = "working directory [required]"),
        make_option(opt_str = c("-r", "--regions"), type = "character", default = NULL, help = "text file of regions to plot, with header, (chr, start, end) [required]"),
        make_option(opt_str = c("-c", "--numCtrl"), type = "integer", default = NULL, help = "total number of control samples [required]"),
        make_option(opt_str = c("-e", "--numExp"), type = "integer", default = NULL, help = "total number of experimental samples [required]"),
        make_option(opt_str = c("-g", "--genome"), type = "character", default = NULL, help = "genome assembly (hg38, hg19, mm10, rn6, rheMac8) [required]"),
        make_option(opt_str = c("-o", "--outprefix"), type = "character", default = NULL, help = "title used in all output files [required]"),
        
        # Optional arguments
        make_option(opt_str = c("--extend"), type = "integer", default = 5000, help = "number of bases to plot on either side of each region [default = 5000]"),
        make_option(opt_str = c("--pctMinCtrl"), type = "double", default = 0.9, help = "minimum percent of control samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--pctMinExp"), type = "double", default = 0.9, help = "minimum percent of experimental samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--mc.cores"), type = "integer", default = 1, help = "cores to use, same as SBATCH -n [default = 1]"),
        make_option(opt_str = c("--estimate.var"), type = "character", default = "same", help = "method to estimate variance for t-test (same, paired, group2) [default = same]"),
        make_option(opt_str = c("--meanDiff_cutoff"), type = "double", default = 0.05, help = "minimum difference between group means for DMRs [default = 0.05]"),
        make_option(opt_str = c("--maxGap"), type = "integer", default = 300, help = "maximum distance between all consecutive CpGs in a DMR [default = 300]"),
        make_option(opt_str = c("--invdensity_cutoff"), type = "integer", default = 300, help = "maximum average distance between consecutive CpGs in a DMR [default = 300]"),
        make_option(opt_str = c("--colorCtrl"), type = "character", default = "3366CC", help = "color for control samples in plots [default = 3366CC]"),
        make_option(opt_str = c("--colorExp"), type = "character", default = "FF3366", help = "color for experimental samples in plots [default = FF3366]")
)
opt_obj <- OptionParser(option_list = option_list, epilogue = "Add DSS_file prefixes after all other arguments, with all control samples first (DSS_files/CTRL01_)", 
                        add_help_option = TRUE, usage = "usage: %prog [arguments]")
opt <- parse_args(object = opt_obj, positional_arguments = c(0, Inf), print_help_and_exit = FALSE)
          
# Test for required arguments
if(opt$options$help){print_help_compact(opt_obj); stop("", call.=FALSE)}
if(is.null(opt$options$chrNum)){print_help_compact(opt_obj); stop("chrNum must be supplied", call.=FALSE)}
if(is.null(opt$options$setwd)){print_help_compact(opt_obj); stop("setwd must be supplied", call.=FALSE)}
if(is.null(opt$options$regions)){print_help_compact(opt_obj); stop("regions must be supplied", call.=FALSE)}
if(is.null(opt$options$numCtrl)){print_help_compact(opt_obj); stop("numCtrl must be supplied", call.=FALSE)}
if(is.null(opt$options$numExp)){print_help_compact(opt_obj); stop("numExp must be supplied", call.=FALSE)}
if(is.null(opt$options$genome)){print_help_compact(opt_obj); stop("genome must be supplied", call.=FALSE)}
if(is.null(opt$options$outprefix)){print_help_compact(opt_obj); stop("outprefix must be supplied", call.=FALSE)}
if(length(opt$args) < 1){print_help_compact(opt_obj); stop("At least 1 DSS file must be supplied", call.=FALSE)}
cat("\n", str(opt))   

# Assign arguments to global variables
chrNum <- as.numeric(opt$options$chrNum)
setwd(as.character(opt$options$setwd))
regions <- read.delim(file = as.character(opt$options$regions), header = FALSE, sep = "\t")
regions <- regions[2:nrow(regions),]
colnames(regions) <- c("chr", "start", "end")
numCtrl <- as.numeric(opt$options$numCtrl)                  
numExp <- as.numeric(opt$options$numExp) 
genome <- as.character(opt$options$genome)
outprefix <- as.character(opt$options$outprefix)
extend <- as.numeric(opt$options$extend)
numMinCtrl <- ceiling(as.numeric(opt$options$pctMinCtrl)*numCtrl)               
numMinExp <- ceiling(as.numeric(opt$options$pctMinExp)*numExp)
mc.cores <- as.numeric(opt$options$mc.cores)                 
estimate.var <- as.character(opt$options$estimate.var)           
meanDiff_cutoff <- as.numeric(opt$options$meanDiff_cutoff)         
maxGap <- as.numeric(opt$options$maxGap)                  
invdensity_cutoff <- as.numeric(opt$options$invdensity_cutoff)       
colorCtrl <- as.character(opt$options$colorCtrl)             
colorCtrl <- paste("#", colorCtrl, sep="")
colorExp <- as.character(opt$options$colorExp)              
colorExp <- paste("#", colorExp, sep="")
DSSprefix <- opt$args   

################################################
# Packages
################################################
cat("[DMRplot] Loading packages\n")
cat("\nlibrary(BiocGenerics)\n"); library(BiocGenerics, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(Biobase)\n"); library(Biobase, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(S4Vectors)\n"); library(S4Vectors, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(matrixStats)\n"); library(matrixStats, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(DelayedArray)\n"); library(DelayedArray, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(bsseq)\n"); library(bsseq, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(DSS)\n"); library(DSS, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(permute)\n"); library(permute, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(GenomicRanges)\n"); library(GenomicRanges, logical.return = TRUE, quietly = TRUE)
cat("\nlibrary(scales)\n"); library(scales, logical.return = TRUE, quietly = TRUE)

#################################################
# DMRplot Pipeline
#################################################

#Set up variables
CTRLgroup <- paste("C",1:numCtrl,sep="")
EXPgroup <- paste("E",1:numExp,sep="")

# Get gene and CpG island bed file names
cat("\n[DMRplot] Loading genes and CpG islands\n")
if(genome == "hg38"){
        CGI_bedfile <- "/share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed"
        Genes_bedfile <- "/share/lasallelab/genomes/hg38/GTF/hg38_RefSeq_Genes.bed"
} else if(genome == "hg19"){
        CGI_bedfile <- "/share/lasallelab/genomes/hg19/GTF/hg19_genome_CGI.bed"
        Genes_bedfile <- "/share/lasallelab/genomes/hg19/GTF/hg19_RefSeq_Genes.bed"
} else if(genome == "mm10"){
        CGI_bedfile <- "/share/lasallelab/genomes/mm10/GTF/mm10_genome_CGI.bed"
        Genes_bedfile <- "/share/lasallelab/genomes/mm10/GTF/mm10_RefSeq_Genes.bed"
} else if(genome == "rn6"){
        CGI_bedfile <- "/share/lasallelab/genomes/rn6/GTF/rn6_genome_CGI.bed"
        Genes_bedfile <- "/share/lasallelab/genomes/rn6/GTF/rn6_RefSeq_Genes.bed"
} else if(genome == "rheMac8"){
        CGI_bedfile <- "/share/lasallelab/genomes/rheMac8/GTF/rheMac8_genome_CGI.bed"                   
        Genes_bedfile <- "/share/lasallelab/genomes/rheMac8/GTF/rheMac8_RefSeq_Genes.bed"
} else{cat(paste("Warning! Gene locations unknown because genome is defined as ",genome,"\n"))}

# Load genes bed file
genome_genes <- read.delim(Genes_bedfile, header=FALSE)
colnames(genome_genes) <- c("chr", "start", "end", "gene_name", "space", "strand")
genome_genes$strand <- as.character(genome_genes$strand)
genome_genes$gene_name <- as.character(genome_genes$gene_name)

# Load CpG islands bed file
genome_CpG_islands <- bed_to_granges(CGI_bedfile)
genome_list <- list(genome_CpG_islands)
names(genome_list) <- c("CGI")

# Load chromosome data
if(genome == "hg38" | genome == "hg19"){chroms = paste("chr",1:22,sep=""); chroms = c(chroms,"chrX","chrY", "chrM")
}else if(genome == "mm10"){chroms = paste("chr",1:19,sep="");chroms = c(chroms,"chrX","chrY", "chrM")  
}else if(genome == "rn6"){chroms = paste("chr",1:20,sep="");chroms = c(chroms,"chrX","chrY", "chrM")
}else if(genome == "rheMac8"){chroms = paste("chr",1:20,sep="");chroms = c(chroms,"chrX","chrY", "chrM")
}else{cat(paste("Warning! Chromosome names unknown because genome is defined as ",genome,"\n"))}

chrom = chroms[chrNum]
cat("\n[DMRplot] Loading", chrom, "DSS files\n")

DSSlist = list(read.table(paste(DSSprefix[1],chrom,".DSS.txt",sep=""), header=TRUE))
for(i in 2:length(DSSprefix))
        DSSlist = c(DSSlist,list(read.table(paste(DSSprefix[i],chrom,".DSS.txt",sep=""), header=TRUE)))

# Make BSobject and smooth
cat("\n[DMRplot] Making BSobject and smoothing\n")
BSobj <- makeBSseqData(DSSlist,c(CTRLgroup,EXPgroup))
rm(DSSlist)
BSobj_smoothed <- BSmooth(BSseq = BSobj, mc.cores = mc.cores, maxGap = 10^8, parallelBy = "sample", ns = 70, 
                         h = 1000, mc.preschedule = FALSE, keep.se = FALSE, verbose = TRUE)    
cat("\n[DMRplot] Completed smoothing\n")

pData <- pData(BSobj_smoothed)
pData$type <- as.character(c(rep(c("Ctrl"), numCtrl), rep(c("Exp") ,numExp)))
pData$col <- c(rep(c(colorCtrl), numCtrl), rep(c(colorExp), numExp))
pData(BSobj_smoothed) <- pData

# Subset BSobject by coverage
cat("\n[DMRplot] Subsetting BSobject by coverage\n")
BSobj_cov <- getCoverage(BSobj_smoothed)
keep_loci <- which(rowSums(BSobj_cov[, BSobj_smoothed$type == "Ctrl"] >= 1) >= numMinCtrl & rowSums(BSobj_cov[, BSobj_smoothed$type == "Exp"] >= 1) >= numMinExp)
BSobj_keep <- BSobj_smoothed[keep_loci,]

# Perform t-tests
cat("\n[DMRplot] Performing t-tests by CpG\n")
BSobj_tstat <- BSmooth.tstat(BSobj_keep, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = estimate.var, 
                                     local.correct = TRUE, qSd = 0.75, k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)    

# Find DMRs
cat("\n[DMRplot] Identifying DMRs\n")
cutoff <- qt(1-0.05/2, numCtrl+numExp-2)
all_DMRs <- dmrFinder(BSobj_tstat, cutoff = c(-cutoff, cutoff), maxGap = maxGap, verbose = FALSE)
if(length(all_DMRs[,1]) > 0){
        all_DMRs$end <- all_DMRs$end + 1 #Add one base to end to include last CpG
        silver_DMRs <- subset(all_DMRs, n >=3 & abs(meanDiff) > meanDiff_cutoff & invdensity <= invdensity_cutoff)
        cat("\n[DMRplot] Found", length(silver_DMRs[,1]), "silver DMRs\n")
}else{cat("Warning, no DMRs found in ", chrom ,"\n")}

# Print output files for silver DMRs
cat("\n[DMRplot] Printing region plots\n")
dmrfilename = paste(outprefix, chrom, "region_plots.pdf", sep = "_")
pdf(file = dmrfilename, width = 6.7, height = 3.25)  #Opens file for figures
plotManyRegions2(BSseq = BSobj_smoothed, regions = regions, extend = extend, addRegions = silver_DMRs,
                lwd = rep(1.1, numCtrl+numExp), verbose = FALSE, BSseqStat = BSobj_tstat, stat.lwd = 1.3,
                stat.ylim = c(-6,6), geneTrack = genome_genes, cex.gene = 1.1, annoTrack = genome_list, cex.anno = 0.8)
dev.off()

cat("\n[DMRplot] Finished!\n\n")


