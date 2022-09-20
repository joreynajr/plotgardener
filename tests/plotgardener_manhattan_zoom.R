library(plotgardener)
library(stringr)
library(plyr)

# Collecting command line arguments
# 1. GWAS summary statistics (TSV format)
# 2. Output path for SVG
#args = commandArgs(trailingOnly=TRUE)
#ss = 'results/main/chiou_2021/processing/finemapping/finemapping.gaulton.tsv'
ss = 'tests/finemapping.gaulton.tsv'
outdir = 'tests/test_finemapping/'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
args = c(ss, outdir)

## Load GWAS data
gwas_ss <- read.table(args[1], sep='\t', header = TRUE)

##### Process LD information
### Plot GWAS data zooming in on chromosome 11
### highlighting a lead SNP, and coloring by LD score
#gwas_ss$LD <- as.numeric(gwas_ss$LD)
### Group LD column into LD ranges
#gwas_ss <- as.data.frame(dplyr::group_by(gwas_ss,
#                                                   LDgrp = cut(
#                                                       gwas_ss$LD,
#                                                       c(0, 0.2, 0.4, 0.6, 0.8, 1))))
#gwas_ss$LDgrp <- addNA(gwas_ss$LDgrp)


plot_gwas_ss <- function(chrom, start, end, ss){

    return
}

colnames(gwas_ss) <- c("chrom", "pos", "snp", "p",
                       "Signal.name", "Major.allele",
                       "Minor.allele", "Unique.ID..hg38.",
                       "Position..hg38.")
gwas_ss["chrom"] <- as.character(gwas_ss[,"chrom"])

# adding chr for visualization
if (sum(c("chr", "Chr") %in% gwas_ss["chrom"]) == 0){
    gwas_ss[,"chrom"] <- paste0("chr",  gwas_ss[,"chrom"])
}


uniq_signals <- unlist(unique(gwas_ss["Signal.name"]))
for (signal_name in uniq_signals){
    
    ## Open an SVG device
    tmp_fn = paste0(signal_name, '.svg')
    tmp_fn = str_replace_all(tmp_fn, ':', '_')
    tmp_fn = file.path(outdir, tmp_fn)
    
    svg(tmp_fn, width = 8.0, height = 3)
    
    ## Create a page
    pageCreate(width = 7.5, height = 2, default.units = "inches")
    
    print(signal_name)
    
    # Dxtract the SNPs from the current signal
    curr_idxs = which(gwas_ss["Signal.name"] == signal_name)
    
    if (length(curr_idxs) == 1){
        next
        
    }
    curr_df = gwas_ss[curr_idxs,]
    example = curr_df[1,]

    ## Extract the lead SNP
    leadSNP_pp <- max(curr_df$p)
    leadSNP <- gwas_ss[which(gwas_ss$p == leadSNP_pp),]$snp
    
    left_pos = min(curr_df$pos)
    left_pos = plyr::round_any(left_pos, 1000, floor) 
    
    right_pos = max(curr_df$pos)
    right_pos = plyr::round_any(right_pos, 1000, ceiling)
    
    chr11_manhattanPlot <- plotManhattanGeneral(
        data = curr_df,
        chrom = example[1,"chrom"],
        chromstart = left_pos,
        chromend = right_pos,
        assembly = "hg19",
        sigLine = FALSE,
        col = "grey",
        lty = 2,
        range = c(0, 1),
        leadSNP = list(
            snp = leadSNP,
            pch = 18,
            cex = 0.75,
            fill = "#7ecdbb",
            fontsize = 8
        ),
        x = 0.5, y = 0, width = 6.5,
        height = 1.5,
        just = c("left", "top"),
        default.units = "inches"
    )
    
    # ## Plot legend for LD scores
    # plotLegend(
    #     legend = c(
    #         "LD Ref Var",
    #         paste("0.4", ">", "r^2", 
    #               "", ">=", "0.2"),
    #         paste("0.2", ">", "r^2", 
    #               "", ">=", "0"),
    #         "no LD data"
    #     ),
    #     fill = c("#7ecdbb", "#37a7db", "#1f4297", "grey"), cex = 0.75,
    #     pch = c(18, 19, 19, 19), border = FALSE, x = 7, y = 0,
    #     width = 1.5, height = 0.6, just = c("right", "top"),
    #     default.units = "inches"
    # )
    
    
    ## Annotate genome label
    annoGenomeLabel(
        plot = chr11_manhattanPlot, x = 0.5, y = 1.5,
        fontsize = 8, scale = "Kb", 
        just = c("left", "top"), default.units = "inches",
        commas = TRUE
    )
    
    ## Annotate y-axis
    annoYaxis(
        plot = chr11_manhattanPlot,
        at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
        axisLine = TRUE, fontsize = 8
    )

    ## Plot y-axis label
    plotText(
        label = "PPA", x = 0.1, y = 0.75, rot = 90,
        fontsize = 8, fontface = "bold", just = "center",
        default.units = "inches"
    )

    ## Hide page guides
    pageGuideHide()
    
    dev.off()
    
    #break
}




