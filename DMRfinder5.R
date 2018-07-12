#!/usr/bin/env Rscript

# Updated 7/12/18
# Charles Mordaunt
# This code is to be used on barbera.genomecenter.ucdavis.edu

cat("\n[DMRfinder]\n\nA pipeline to identify differentially-methylated regions from whole-genome bisulfite sequencing data\n\n")

#################################################
# Functions 
#################################################

cat("\n[DMRfinder] Loading functions\n\n")

write.bed = function(df,file,header="",name="name") {
# Written by Keith Dunaway
# Generalized write bed function
# Example Output:
#  chr1    10469   10470   0.00-1  0       +       0       0       0,0,0
# Area stat as name
# color for direction (hyper = red, hypo = blue)
        colnum = 3
        if(name %in% colnames(df)) {namedat = df[,name];colnum = 4} else{namedat = rep(".",length(df[,1]))}
        if("score" %in% colnames(df)) {score = df$score;colnum = 5} else{score = rep("0",length(df[,1]))}
        if("strand" %in% colnames(df)) {strand = df$strand;colnum = 6} else{strand = rep("+",length(df[,1]))}
        if("thickStart" %in% colnames(df)) {thickStart = df$thickStart;colnum = 7} else{thickStart = rep("0",length(df[,1]))}
        if("thickEnd" %in% colnames(df)) {thickEnd = df$thickEnd;colnum = 8} else{thickEnd = rep("0",length(df[,1]))}
        if("itemRgb" %in% colnames(df)) {itemRgb = df$itemRgb;colnum = 9} else{itemRgb = rep("0,0,0",length(df[,1]))}
        if("blockCount" %in% colnames(df)) {blockCount = df$blockCount;colnum = 10} else{blockCount = rep("0",length(df[,1]))}
        if("blockSizes" %in% colnames(df)) {blockSizes = df$blockSizes;colnum = 11} else{blockSizes = rep("0",length(df[,1]))}
        if("blockStarts" %in% colnames(df)) {blockStarts = df$blockStarts;colnum = 12} else{blockStarts = rep("0",length(df[,1]))}
        subdf = cbind.data.frame(df$chr,df$start,df$end,namedat,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)
        if(header == ""){
                if(length(subdf[,1]) == 1){
                        write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE, eol="\t")
                }else{
                        write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE, eol="\n")
                }
        }else{
                if(length(subdf[,1]) == 1){
                        write(header,file=file)
                        write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE,append=TRUE, eol="\t")
                }else{
                        write(header,file=file)
                        write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE,append=TRUE, eol="\n")
                }
        }
}

write.dmrs_bed = function(df,file,trackname,genome) {
# Written by Keith Dunaway
# Wrapper function for writing DMR info data.frames as bed files, Area stat as name, color for direction (hyper = red, hypo = blue)
# Takes in DMR info data.frame and formats it for input into write.bed function        
        df$itemRgb = ifelse(df$direction == "hypo","0,0,255","255,0,0")
        headerline = paste("track name=",trackname," description=",trackname," useScore=0 itemRgb=On genome=",genome,sep="")
        write.bed(df,file,header=headerline,name="areaStat")
}

print_help_compact <- function (object) {
# Originally written by Trevor Davis as part of optparse package
# Modified by Charles Mordaunt
# Prints help for an Option Parser object, more compact than print_help()
        cat(object@usage, fill = TRUE)
        cat(object@description, fill = TRUE)
        cat("Required arguments:", sep = "\n")
        options_list <- object@options
        for (ii in 1:6) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat("\nOptional arguments:", sep = "\n")
        for (ii in 7:16) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat("\nOutput arguments:", sep = "\n")
        for (ii in 17:26) {
                option <- options_list[[ii]]
                cat("\t")
                if (!is.na(option@short_flag)) cat(option@short_flag, ", ", sep = "")
                if (!is.null(option@long_flag)) cat(option@long_flag, " = ", option@help, sep = "")
                cat("\n")
        }
        cat(object@epilogue, fill = TRUE, sep = "")
        return(invisible(NULL))
}

permuteAll <- function(nperm, design) {
# Written by Kasper Daniel Hansen        
# Downloaded from https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R
# Creates an matrix of permutations based on design (number of samples), and nperm (number of permutations)
        message(sprintf("[permuteAll] performing %d unrestricted permutations of the design matrix\n", nperm))
        CTRL <- how(nperm = nperm)
        idxMatrix <- shuffleSet(n = design, control = CTRL)
}

subsetByMatrix <- function(vec, mat) {
# Written by Kasper Daniel Hansen        
# Downloaded from https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R
# Assigns sample names to numbers in permutations
        apply(mat, 2, function(xx) vec[xx])
}

getNullDistribution_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var, local.correct = TRUE, cutoff, stat = "tstat.corrected", maxGap, mc.cores) {
# Originally written by Kasper Daniel Hansen        
# Downloaded from https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R
# Modified by Charles Mordaunt
# For each permutation, creates tstat object and finds DMRs, returns a list
        stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
        message(sprintf("[getNullDistribution_BSmooth.tstat] performing %d permutations", nrow(idxMatrix1)-1)) 
        nullDist <- lapply(1:nrow(idxMatrix1), function(ii) {
                ptime1 <- proc.time()
                BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var, group1 = idxMatrix1[ii,], group2 = idxMatrix2[ii,], local.correct = local.correct, 
                                          maxGap = 10^8, verbose = FALSE, mc.cores = mc.cores, qSd = 0.75, k = 101)
                dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = maxGap, verbose = FALSE)
                dmrs0$end <- dmrs0$end + 1 #Add 1 base to end to include last CpG
                ptime2 <- proc.time()
                stime <- (ptime2 - ptime1)[3]
                message(sprintf("[getNullDistribution_BSmooth.tstat] completing permutation %d in %.1f sec", ii-1, stime))
                dmrs0
        })
        nullDist
}

subsetDmrs <- function(xx, meanDiff, invdensity) {
# Originally written by Kasper Daniel Hansen        
# Downloaded from https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R   
# Modified by Charles Mordaunt
# Subset Null DMRs for CpGs, meanDiff, and invdensity
        if(is.null(xx) || is(xx, "try-error")) return(NULL)
        if(length(xx$end) < 1) return(NULL)
        out <- xx[ xx[,"n"] >= 3 & abs(xx[, "meanDiff"]) > meanDiff & xx[, "invdensity"] <= invdensity, ]
        if(nrow(out) == 0) return(NULL)
        out
}

getFWER <- function(null) {
# Originally written by Kasper Daniel Hansen        
# Downloaded from https://github.com/kasperdanielhansen/bsseq/blob/master/R/permutations.R
# Modified by Charles Mordaunt
# Calculates FWER from subsetted null DMR list
        reference <- null[[1]]
        null <- null[-1]
        null <- null[!sapply(null, is.null)]
        if(!is.null(unlist(null))){
                better <- sapply(1:nrow(reference), function(ii) {
                        areaStat <- abs(reference$areaStat[ii])
                        n <- reference$n[ii]
                        out <- sapply(null, function(nulldist) {
                                any(abs(nulldist$areaStat) >= areaStat & nulldist$n >= n)
                        })
                        sum(out)
                })
                return(better)
        } else {
                return(rep(0, length(reference[,1])))
        }
}

plotGeneTrack2 <- function (gr, geneTrack, cex) {
# Originally written by Kasper Daniel Hansen as part of the bsseq package
# Modified by Charles Mordaunt
# Takes data.frame of gene locations and plots them along with DMRs
        geneTrack_gr <- makeGRangesFromDataFrame(geneTrack)
        ol <- findOverlaps(geneTrack_gr, gr)
        genes <- geneTrack[queryHits(ol), ]
        plot(start(gr), rep(1, length(start(gr))), type = "n", xaxt = "n", yaxt = "n", bty = "n", ylim = c(-1.5, 1.5), xlim = c(start(gr), end(gr)), 
             xlab = "", ylab = "", cex.lab = 4, lheight = 2, cex.axis = 1)
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
                        } else {
                                apply(Exons, 1, function(x) polygon(c(x[1], x[2], x[2], x[1]), c(-0.25, -0.25, -0.65, -0.65), col = "#000099"))
                                text(x = (max(start(gr), min(ES)) + min(end(gr), max(EE)))/2, y = -1, labels = genes$gene_name[g], cex = cex, font = 3)
                        }
                }
        }
}

plotRegion2 <- function (BSseq, region = NULL, extend = 0, main = "", addRegions = NULL, annoTrack = NULL, cex.anno = 1, geneTrack = NULL, cex.gene = 1.5, 
                         col = NULL, lty = NULL, lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected", stat.col = "black", stat.lwd = 1, stat.lty = 1, 
                         stat.ylim = c(-8, 8), mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, 
                         highlightMain = FALSE) {
# Originally written by Kasper Daniel Hansen as part of the bsseq package
# Modified by Charles Mordaunt
# Plots smoothed methylation by sample and t-statistic, along with DMR and gene locations, for a single region
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
        .plotSmoothData2(BSseq = BSseq, region = region, extend = extend, addRegions = addRegions, col = col, lty = lty, lwd = lwd, regionCol = regionCol, 
                         addTicks = addTicks, addPoints = addPoints, pointsMinCov = pointsMinCov, highlightMain = highlightMain)
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
                plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n", ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "", 
                     ylab = stat.type, cex.lab = 1.3)
                axis(side = 2, at = c(-5, 0, 5), cex.axis = 1.2)
                abline(h = 0, col = "grey60")
                bsseq:::.bsPlotLines(start(BSseqStat), stat.values, lty = stat.lty, col = stat.col, lwd = stat.lwd, plotRange = c(start(gr), end(gr)))
                bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(-5.75, 5.75), regionCol = regionCol, highlightMain = highlightMain)
        }
        if (!is.null(annoTrack)) plotAnnoTrack2(gr, annoTrack, cex.anno)
        if (!is.null(geneTrack)) plotGeneTrack2(gr, geneTrack, cex.gene)
        if (!is.null(main)) {
                main <- bsseq:::.bsPlotTitle(gr = region, extend = extend, main = main, mainWithWidth = mainWithWidth)
                mtext(side = 3, text = main, outer = TRUE, cex = 0.9, padj = 0)
        }
        return(invisible(NULL))
}

plotManyRegions2 <- function (BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL, annoTrack = NULL, cex.anno = 1, geneTrack = NULL, cex.gene = 1.5, 
                              col = NULL, lty = NULL, lwd = NULL, BSseqStat = NULL, stat = "tstat.corrected", stat.col = "black", stat.lwd = 1, stat.lty = 1,
                              stat.ylim = c(-8, 8), mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE, pointsMinCov = 5, 
                              highlightMain = FALSE, verbose = TRUE) {
# Originally written by Kasper Daniel Hansen as part of the bsseq package
# Modified by Charles Mordaunt
# Wrapper function for plotRegion2, which takes data.frame or GRanges object and plots every region.
        cat("[plotManyRegions] preprocessing ...")
        if (!is.null(regions)) {
                if (is(regions, "data.frame")) gr <- data.frame2GRanges(regions, keepColumns = FALSE)
                else gr <- regions
                if (!is(gr, "GRanges")) stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
        } else gr <- granges(BSseq)
        gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
        BSseq <- subsetByOverlaps(BSseq, gr)
        if (length(start(BSseq)) == 0) stop("No overlap between BSseq data and regions")
        if (!is.null(main) && length(main) != length(gr)) main <- rep(main, length = length(gr))
        cat("done\n")
        for (ii in seq(along = gr)) {
                if (verbose) cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n", ii, nrow(regions)))
                plotRegion2(BSseq = BSseq, region = regions[ii, ], extend = extend, col = col, lty = lty, lwd = lwd, main = main[ii], 
                            BSseqStat = BSseqStat, stat = stat, stat.col = stat.col, stat.lwd = stat.lwd, stat.lty = stat.lty, stat.ylim = stat.ylim, 
                            addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth, annoTrack = annoTrack, 
                            cex.anno = cex.anno, geneTrack = geneTrack, cex.gene = cex.gene, addTicks = addTicks, addPoints = addPoints, 
                            pointsMinCov = pointsMinCov, highlightMain = highlightMain)
        }
}

.plotSmoothData2 <- function (BSseq, region, extend, addRegions, col, lty, lwd, regionCol, addTicks, addPoints, pointsMinCov, highlightMain) {
# Originally written by Kasper Daniel Hansen as part of the bsseq package
# Modified by Charles Mordaunt
# Plots smoothed methylation data for plotRegion2
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
        plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n", ylim = c(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.05), xlim = c(start(gr), end(gr)), 
             xlab = "", ylab = "Methylation", cex.lab = 1.3)
        axis(side = 2, at = c(ymin, ymid, ymax), cex.axis = 1.2, labels = c(ymin*100, ymid*100, ymax*100))
        if (addTicks) rug(positions, ticksize = 0.08)
        bsseq:::.bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(ymin, ymax), regionCol = regionCol, highlightMain = highlightMain)
        if (addPoints) {
                sapply(1:ncol(BSseq), function(sampIdx) {
                        abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
                })
        }
        sapply(1:ncol(BSseq), function(sampIdx) {
                bsseq:::.bsPlotLines(positions, smoothPs[, sampIdx], col = alpha(colEtc$col[sampIdx], 0.7), lty = colEtc$lty[sampIdx], 
                                     lwd = colEtc$lwd[sampIdx], plotRange = c(start(gr), end(gr)))
        })
        if (addPoints) {
                sapply(1:ncol(BSseq), function(sampIdx) {
                        bsseq:::.bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx], col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov)
                })
        }
}

plotAnnoTrack2 <- function (gr, annoTrack, cex) {
# Originally written by Kasper Daniel Hansen as part of the bsseq package
# Modified by Charles Mordaunt
# Takes GRanges object of genomic annotations and plots it as part of plotRegion2
        if (!all(sapply(annoTrack, function(xx) is(xx, "GRanges")))) stop("all elements in 'annoTrack' needs to be 'GRanges'")
        plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n", ylim = c(0.5, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), 
             xlab = "", ylab = "")
        lapply(seq(along = annoTrack), function(ii) {
                jj <- length(annoTrack) + 1 - ii
                ir <- subsetByOverlaps(annoTrack[[ii]], gr)
                if (length(ir) > 0) rect(start(ir) - 0.5, jj - 0.3, end(ir), jj + 0.3, col = "#006600", lwd = 0.8)
                mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1, cex = cex, adj = 0.6)
        })
}

bed_to_granges <- function(file){
# Written by Dave Tang
# Downloaded from https://github.com/davetang/bedr/blob/master/R/bed_to_granges.R
# Converts bed file to GRanges object for plotAnnoTrack2
        df <- read.table(file, header=F, stringsAsFactors=F)
        if(length(df) > 6) df <- df[,-c(7:length(df))]
        if(length(df)<3) stop("File has less than 3 columns")
        header <- c('chr','start','end','id','score','strand')
        names(df) <- header[1:length(names(df))]
        if('strand' %in% colnames(df)) df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
        if(length(df)==3) gr <- with(df, GRanges(chr, IRanges(start, end)))
        else if (length(df)==4) gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
        else if (length(df)==5) gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
        else if (length(df)==6) gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
        return(gr)
}

#################################################
# Global Variables
#################################################

# Get arguments from bash script
library(optparse)
cat("[DMRfinder] Getting arguments from bash script\n\n")
option_list <- list(
        # Required arguments
        make_option(opt_str = c("-n", "--chrNum"), type = "integer", default = NULL, help = "chromosome number"),
        make_option(opt_str = c("-d", "--setwd"), type = "character", default = NULL, help = "working directory"),
        make_option(opt_str = c("-c", "--numCtrl"), type = "integer", default = NULL, help = "number of control samples"),
        make_option(opt_str = c("-e", "--numExp"), type = "integer", default = NULL, help = "number of experimental samples"),
        make_option(opt_str = c("-a", "--genome"), type = "character", default = NULL, help = "genome assembly (hg38, hg19, mm10, rn6, rheMac8)"),
        make_option(opt_str = c("-o", "--outprefix"), type = "character", default = NULL, help = "title used in output files"),
        
        # Optional arguments
        make_option(opt_str = c("--pctMinCtrl"), type = "double", default = 0.9, help = "minimum percent of control samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--pctMinExp"), type = "double", default = 0.9, help = "minimum percent of experimental samples with 1 read at CpG [default = 0.9]"),
        make_option(opt_str = c("--mc.cores"), type = "integer", default = 1, help = "cores to use, same as SBATCH -n [default = 1]"),
        make_option(opt_str = c("--estimate.var"), type = "character", default = "same", help = "method to estimate variance for t-test (same, paired, group2) [default = same]"),
        make_option(opt_str = c("--meanDiff_cutoff"), type = "double", default = 0.05, help = "minimum difference between group means for DMRs [default = 0.05]"),
        make_option(opt_str = c("--maxGap"), type = "integer", default = 300, help = "maximum distance between all consecutive CpGs in a DMR [default = 300]"),
        make_option(opt_str = c("--invdensity_cutoff"), type = "integer", default = 300, help = "maximum average distance between consecutive CpGs in a DMR [default = 300]"),
        make_option(opt_str = c("--colorCtrl"), type = "character", default = "3366CC", help = "color for control samples in plots [default = 3366CC]"),
        make_option(opt_str = c("--colorExp"), type = "character", default = "FF3366", help = "color for experimental samples in plots [default = FF3366]"),
        make_option(opt_str = c("--nperm"), type = "integer", default = 1000, help = "number of permutations to perform for FWER estimation [default = 1000]"),
        
        # Output arguments
        make_option(opt_str = c("--gold_bed"), type = "logical", default = TRUE, help = "Output bed file of locations for gold DMRs [default = TRUE]"),
        make_option(opt_str = c("--silver_bed"), type = "logical", default = TRUE, help = "Output bed file of locations for silver DMRs [default = TRUE]"),
        make_option(opt_str = c("--silver_info"), type = "logical", default = TRUE, help = "Output txt file of info for silver DMRs [default = TRUE]"),
        make_option(opt_str = c("--gold_plots"), type = "logical", default = TRUE, help = "Output pdf file of plots for gold DMRs [default = TRUE]"),
        make_option(opt_str = c("--silver_plots"), type = "logical", default = FALSE, help = "Output pdf file of plots for silver DMRs [default = FALSE]"),
        make_option(opt_str = c("--background"), type = "logical", default = FALSE, help = "Output bed file of background DMRs [default = FALSE]"),
        make_option(opt_str = c("--meth"), type = "logical", default = FALSE, help = "Output table of smoothed methylation for each sample at each silver DMR [default = FALSE]"),
        make_option(opt_str = c("--cov"), type = "logical", default = FALSE, help = "Output table of total coverage for each sample at each DMR [default = FALSE]"),
        make_option(opt_str = c("--background_cov"), type = "logical", default = FALSE, help = "Output table of total coverage for each sample at each background DMR [default = FALSE]"),
        make_option(opt_str = c("--CpGs"), type = "logical", default = FALSE, help = "Output bed file of tested CpGs [default = FALSE]")
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
genome <- as.character(opt$options$genome)
outprefix <- as.character(opt$options$outprefix)
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
nperm <- as.numeric(opt$options$nperm) 
gold_bed <- as.logical(opt$options$gold_bed)
silver_bed <- as.logical(opt$options$silver_bed)
silver_info <- as.logical(opt$options$silver_info)
gold_plots <- as.logical(opt$options$gold_plots)
silver_plots <- as.logical(opt$options$silver_plots)
background <- as.logical(opt$options$background)
meth <- as.logical(opt$options$meth)
cov <- as.logical(opt$options$cov)
background_cov <- as.logical(opt$options$background_cov)
CpGs <- as.logical(opt$options$CpGs)
DSSprefix <- opt$args   

################################################
# Packages
################################################
cat("[DMRfinder] Loading packages\n")
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
# DMRfinder Pipeline Setup
#################################################

#Set up group variables
CTRLgroup <- paste("C",1:numCtrl,sep="")
EXPgroup <- paste("E",1:numExp,sep="")
groups <- c(EXPgroup, CTRLgroup)

# Create Combination Matrices
# Combinations are created differently based on number possible
#       If < nperm possible (usually 1000), all combinations are created
#       If > nperm, but < 1e8, nperm possibilities are sampled from all combinations
#       If > 1e8, nperm permutations are performed (identical combinations are possible, but unlikely)
cat("\n[DMRfinder] Creating combination matrices\n")
set.seed(1)
if(choose(numCtrl+numExp, numExp) < 1e8){
        idxMatrix <- NULL
        temp <- combn(groups, length(EXPgroup), FUN = NULL, simplify = TRUE)
        temp <- temp[,2:length(temp[1,])]
        if(numCtrl == numExp) temp <- temp[,1:length(temp[1,])-1] # Remove inverse of reference combination
        if(choose(numCtrl+numExp, numExp) > nperm) {
                sample <- sample(x = 1:length(temp[1,]), size = nperm, replace = FALSE)
                temp <- temp[,sample]
                rm(sample)
        }
        for(i in 1:length(temp[1,])){
                t <- c(temp[,i],setdiff(groups,temp[,i]))
                idxMatrix <- rbind(idxMatrix,t)
        }
        idxMatrix1 <- rbind(EXPgroup, idxMatrix[,1:numExp])
        idxMatrix2 <- rbind(CTRLgroup, idxMatrix[,(numExp+1):(numExp+numCtrl)])
        rm(idxMatrix, temp, i, t)
} else {
        idxMatrix <- permuteAll(nperm = nperm, design = numCtrl+numExp)
        idxMatrixsub <- subsetByMatrix(c(CTRLgroup, EXPgroup), idxMatrix)
        idxMatrix1 <- rbind(EXPgroup, idxMatrixsub[,1:numExp])
        idxMatrix2 <- rbind(CTRLgroup, idxMatrixsub[,(numExp+1):(numExp+numCtrl)])
        rm(idxMatrix,idxMatrixsub)
}

# Get gene and CpG island bed file names
cat("\n[DMRfinder] Loading genes and CpG islands\n")
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

#################################################
# Main DMRfinder Pipeline
#################################################

# Load chromosome DSS files
chrom = chroms[chrNum]
cat("\n[DMRfinder] Loading", chrom, "DSS files\n")
DSSlist = list(read.table(paste(DSSprefix[1],chrom,".DSS.txt",sep=""), header=TRUE))
for(i in 2:length(DSSprefix)) DSSlist = c(DSSlist,list(read.table(paste(DSSprefix[i],chrom,".DSS.txt",sep=""), header=TRUE)))

# Make BSobject and smooth
cat("\n[DMRfinder] Making BSobject and smoothing\n")
BSobj <- makeBSseqData(DSSlist,c(CTRLgroup,EXPgroup))
rm(DSSlist)
BSobj_smoothed <- BSmooth(BSseq = BSobj, mc.cores = mc.cores, maxGap = 10^8, parallelBy = "sample", ns = 70, h = 1000, mc.preschedule = FALSE, 
                          keep.se = FALSE, verbose = TRUE)    
cat("\n[DMRfinder] Completed smoothing\n")

# Add group assignments and color to smoothed BSobject
pData <- pData(BSobj_smoothed)
pData$type <- as.character(c(rep(c("Ctrl"), numCtrl), rep(c("Exp") ,numExp)))
pData$col <- c(rep(c(colorCtrl), numCtrl), rep(c(colorExp), numExp))
pData(BSobj_smoothed) <- pData

# Subset BSobject by coverage
cat("\n[DMRfinder] Subsetting BSobject by coverage\n")
BSobj_cov <- getCoverage(BSobj_smoothed)
keep_loci <- which(rowSums(BSobj_cov[, BSobj_smoothed$type == "Ctrl"] >= 1) >= numMinCtrl & rowSums(BSobj_cov[, BSobj_smoothed$type == "Exp"] >= 1) >= numMinExp)
BSobj_keep <- BSobj_smoothed[keep_loci,]

# Perform t-tests
cat("\n[DMRfinder] Performing t-tests by CpG\n")
BSobj_tstat <- BSmooth.tstat(BSobj_keep, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = estimate.var, local.correct = TRUE, qSd = 0.75, k = 101, 
                             maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)    

# Output tested CpGs in bed file
if(CpGs){
        CpG_keep <- data.frame(chr = rep(as.character(BSobj_keep@rowRanges@seqnames@values), BSobj_keep@rowRanges@seqnames@lengths), 
                               start = as.numeric(BSobj_keep@rowRanges@ranges@start),
                               end = as.numeric(BSobj_keep@rowRanges@ranges@start)+1)
        cat("\n[DMRfinder] Printing locations for", length(CpG_keep[,1]), "tested CpGs\n")
        write.table(CpG_keep, paste(outprefix, chrom, "CpGs.bed", sep="_"), sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Find DMRs
cat("\n[DMRfinder] Identifying DMRs\n")
cutoff <- qt(1-0.05/2, numCtrl+numExp-2)
all_DMRs <- dmrFinder(BSobj_tstat, cutoff = c(-cutoff, cutoff), maxGap = maxGap, verbose = FALSE)
if(length(all_DMRs[,1]) > 0){
        all_DMRs$end <- all_DMRs$end + 1 #Add one base to end to include last CpG
        silver_DMRs <- subset(all_DMRs, n >=3 & abs(meanDiff) > meanDiff_cutoff & invdensity <= invdensity_cutoff)
        cat("\n[DMRfinder] Found", length(silver_DMRs[,1]), "silver DMRs\n")
        
        # Do permutation testing to estimate FWER
        if(nperm > 0){
                if(length(silver_DMRs[,1]) > 0){
                        cat("\n[DMRfinder] Performing permutation testing to estimate FWER\n")
                        null_DMRs <- getNullDistribution_BSmooth.tstat(BSobj_keep, idxMatrix1, idxMatrix2, cutoff = c(-cutoff, cutoff), 
                                                                       mc.cores = mc.cores, maxGap = maxGap, estimate.var = estimate.var, local.correct = TRUE)
                        null_silver_DMRs <- lapply(null_DMRs, subsetDmrs, meanDiff = meanDiff_cutoff, invdensity = invdensity_cutoff)
                        silver_DMRs$FWER <- getFWER(null_silver_DMRs)
                        if(max(silver_DMRs$FWER) == 0){
                                silver_DMRs$Rel_FWER <- rep(0, length(silver_DMRs[,1]))  
                        } else {
                                silver_DMRs$Rel_FWER <- silver_DMRs$FWER / max(silver_DMRs$FWER)           
                        }                
                        gold_DMRs <- subset(silver_DMRs, Rel_FWER < 0.05)
                        if(length(gold_DMRs[,1]) > 0){
                                cat("\n[DMRfinder] Found", length(gold_DMRs[,1]), "gold DMR(s)!\n")
                                if(gold_bed){
                                        cat("\n[DMRfinder] Printing gold DMR locations\n")
                                        write.dmrs_bed(gold_DMRs, paste(outprefix, chrom, "gold_DMRs.bed", sep = "_"), paste(outprefix, chrom, "gold_DMRs", sep = "_"), genome)
                                }
                                if(gold_plots){
                                        cat("\n[DMRfinder] Printing gold DMR plots\n")
                                        dmrfilename = paste(outprefix, chrom,"gold_DMR_plots.pdf",sep="_")
                                        pdf(file = dmrfilename, width = 6.7, height = 3.25)
                                        plotManyRegions2(BSseq = BSobj_smoothed, regions = gold_DMRs, extend = 5000, addRegions = silver_DMRs,
                                                         lwd = rep(1.1, numCtrl+numExp), verbose = FALSE, BSseqStat = BSobj_tstat, stat.lwd = 1.3,
                                                         stat.ylim = c(-6,6), geneTrack = genome_genes, cex.gene = 1.1, annoTrack = genome_list, cex.anno = 0.8)   
                                        dev.off()
                                }
                        } else{cat("\n[DMRfinder] Warning, no gold DMRs found in", chrom ,"\n")}
                }else{cat("\n[DMRfinder] Warning, no silver DMRs found in", chrom ,"\n")}
        }
        
        # Print output files for silver DMRs
        if(length(silver_DMRs[,1]) > 0){
                if(silver_bed){
                        cat("\n[DMRfinder] Printing silver DMR locations\n")
                        write.dmrs_bed(silver_DMRs, paste(outprefix, chrom, "silver_DMRs.bed", sep = "_"), paste(outprefix, chrom, "silver_DMRs", sep="_"), genome)
                }
                if(silver_plots){
                        cat("\n[DMRfinder] Printing silver DMR plots\n")
                        dmrfilename = paste(outprefix, chrom, "silver_DMR_plots.pdf", sep = "_")
                        pdf(file = dmrfilename, width = 6.7, height = 3.25)  #Opens file for figures
                        plotManyRegions2(BSseq = BSobj_smoothed, regions = silver_DMRs, extend = 5000, addRegions = silver_DMRs,
                                         lwd = rep(1.1, numCtrl+numExp), verbose = FALSE, BSseqStat = BSobj_tstat, stat.lwd = 1.3,
                                         stat.ylim = c(-6,6), geneTrack = genome_genes, cex.gene = 1.1, annoTrack = genome_list, cex.anno = 0.8)
                        dev.off()
                }
                if(silver_info){
                        cat("\n[DMRfinder] Printing silver DMR info\n")
                        silver_DMR_info <- silver_DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", 
                                                          "group2.mean", "group1.mean", "meanDiff", "direction", "FWER", "Rel_FWER")]
                        colnames(silver_DMR_info) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", 
                                                       "Ctrl_mean", "Exp_mean", "meanDiff", "direction", "FWER", "Rel_FWER")
                        write.table(silver_DMR_info, paste(outprefix, chrom, "silver_DMR_info.txt", sep = "_"), sep  = "\t")
                }
        }
        
        # Output smoothed methylation in each DMR for each sample
        if(length(silver_DMRs[,1]) > 0){
                if(meth){
                        cat("\n[DMRfinder] Printing silver DMR smoothed methylation\n")
                        DMR_meth <- data.frame(getMeth(BSseq = BSobj_keep, regions = silver_DMRs, type = "smooth", what = "perRegion"))
                        DMR_meth <- round(DMR_meth, 4)
                        sample_names <- gsub(".*/", "", DSSprefix)
                        sample_names <- gsub("_", "", sample_names)
                        colnames(DMR_meth) <- sample_names
                        DMR_meth <- cbind(silver_DMRs[,c("chr", "start", "end")], DMR_meth)
                        write.table(DMR_meth, paste(outprefix, chrom,"silver_DMR_methylation.txt", sep="_"), sep ="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                }
        
        # Output coverage in each DMR for each sample
                if(cov){
                        cat("\n[DMRfinder] Printing silver DMR coverage\n")
                        DMR_cov <- data.frame(getCoverage(BSseq = BSobj_keep, regions = silver_DMRs, type = "Cov", what = "perRegionTotal"))
                        sample_names <- gsub(".*/", "", DSSprefix)
                        sample_names <- gsub("_", "", sample_names)
                        colnames(DMR_cov) <- sample_names
                        DMR_cov <- cbind(silver_DMRs[,c("chr", "start", "end")], DMR_cov)
                        write.table(DMR_cov, paste(outprefix, chrom, "silver_DMR_coverage.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                }
        }
        
}else{cat("Warning, no DMRs found in ", chrom ,"\n")}

# Make Background DMR bed file
if(background){
        cat("\n[DMRfinder] Identifying background DMRs\n")
        all_background_DMRs <- dmrFinder(BSobj_tstat, cutoff = c(-0, 0), maxGap = maxGap, verbose = FALSE)  
        if(length(all_background_DMRs[,1]) > 0){
                all_background_DMRs$end <- all_background_DMRs$end + 1 #Add one base to end to include last CpG
                silver_background_DMRs <- subset(all_background_DMRs, n >=3 & invdensity <= invdensity_cutoff)
                cat("\n[DMRfinder] Found", length(silver_background_DMRs[,1]), "background DMRs\n")
        }else cat("Warning, no background DMRs found in ", chrom ,"\n")
        if(length(silver_background_DMRs[,1]) > 0){
                cat("\n[DMRfinder] Printing background DMR locations\n")
                headerline <- paste("track name=", paste(outprefix, chrom, "background_DMRs",sep = "_")," description=", paste(outprefix, chrom, "background_DMRs", sep = "_"),
                                    " useScore=0 itemRgb=Off genome=", genome, sep = "")
                write(headerline, paste(outprefix, chrom, "background_DMRs.bed", sep = "_"))
                write.table(silver_background_DMRs[,c("chr", "start", "end")], paste(outprefix, chrom, "background_DMRs.bed", sep = "_"), append = TRUE, quote = FALSE, sep = "\t",
                            row.names = FALSE, col.names = FALSE)
                if(background_cov){
                        cat("\n[DMRfinder] Printing background DMR coverage\n")
                        back_cov <- data.frame(getCoverage(BSseq = BSobj_keep, regions = silver_background_DMRs, type = "Cov", what = "perRegionTotal"))
                        sample_names <- gsub(".*/", "", DSSprefix)
                        sample_names <- gsub("_", "", sample_names)
                        if(length(silver_background_DMRs[,1]) == 1){
                                back_cov <- t(as.matrix(back_cov))
                                back_cov <- cbind(silver_background_DMRs[,c("chr", "start", "end")], back_cov)
                                colnames(back_cov)[4:length(back_cov[1,])] <- sample_names
                                write.table(back_cov, paste(outprefix, chrom, "background_DMR_coverage.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                        } else {
                                colnames(back_cov) <- sample_names
                                back_cov <- cbind(silver_background_DMRs[,c("chr", "start", "end")], back_cov)
                                write.table(back_cov, paste(outprefix, chrom, "background_DMR_coverage.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
                        }
                }
        }else cat("Warning, no silver background DMRs found in ", chrom ,"\n")
}

# Cleanup memory
rm(list=ls())
cat("\n[DMRfinder] Finished!\n\n")


