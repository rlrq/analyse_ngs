#!/usr/bin/Rscript

library(optparse)

option_list <- list(
    make_option(c("-i", "--in"), type = "character", dest = "fin",
                help = "file output by intersect step of analyse_ngs (required)", metavar = "FILE"),
    make_option(c("-o", "--out"), type = "character", dest = "fout",
                help = "path to output file (required)", metavar = "FILE"),
    make_option(c("--collapse"), action = "store_true", dest = "collapse", default = FALSE,
                help = "raise this to collapse 'gene' in to a single row per query sequence")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

## check arguments
if (is.na(args$fin)){
    stop("Input file must be provided. Use '--in <PATH>'.")
}
if (is.na(args$fou)){
    stop("Ouput file must be provided. Use '--out <PATH>'")
}

library(magrittr)

## get best gene(s) by bitscore per derep seq
df.intersect <- read.table(args$fin, header = TRUE, sep = '\t')
df.best <- df.intersect %>%
    dplyr::group_by(qseqid) %>%
    dplyr::slice_max(bitscore, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, bitscore, gene)

## collapse gene column
if (args$collapse){
    df.best <- df.best %>%
        dplyr::group_by(qseqid, bitscore) %>%
        dplyr::summarise(gene = paste0(gene, collapse = ',')) %>%
        dplyr::ungroup()
}

## write
write.table(df.best, file = args$fout, quote = FALSE, sep = '\t', row.names = FALSE, col.name = TRUE)
