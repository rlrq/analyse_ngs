#!/usr/bin/Rscript

library(optparse)

## pipeline adapted from: https://astrobiomike.github.io/amplicon/dada2_workflow_ex
## (also, zzOtherzz/ChenQi_16S/src/analysis_16S.R)
## (also, multi_ko/src/derep_fastq.R)
## (also, zzOtherzz/YiYun/ploop_primer/src/derep_fastq.R)

option_list <- list(
    make_option(c("-f", "--fastq"), type = "character", dest = "dir_fastq",
                help = "path to directory of fastq files (required)", metavar = "DIR"),
    make_option(c("-o", "--out"), type = "character", dest = "dir_output",
                help = "path to output directory (required)", metavar = "DIR"),
    make_option(c("--gz"), action = "store_true", dest = "gz", default = FALSE,
                help = "raise this flag if fastq files are gzipped"),
    make_option(c("-p", "--prefix"), type = "character", dest = "prefix", default = "derepTrim",
                help = "prefix of output files/directories", metavar = "STRING"),
    make_option(c("-t", "--truncation"), type = "integer", dest = "trunc_len", default = 70,
                help = "trim reads to first X bp", metavar = "INTEGER"),
    make_option(c("-s", "--fastq-suffix"), type = "character", dest = "fastq_suffix_pattern",
                default = "_L001_R\\d_001",
                help = paste0("fastq file suffix regex;",
                              " use 'R\\d' pattern for R1/R2 matching;",
                              " do not include extensions '.fastq' and/or '.gz'")),
    ## dada2 function options
    make_option(c("--maxEE"), type = "numeric", dest = "maxEE", default = 2,
                help = "maxEE parameter for dada2::filterAndTrim"),
    make_option(c("--rm-phix"), type = "logical", dest = "rm_phix", default = TRUE,
                help = "rm.phix parameter for dada2::filterAndTrim"),
    make_option(c("-m", "--multithread"), type = "logical", dest = "multithread", default = TRUE,
                help = "multithread parameter for dada2::learnErrors")
)
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

## check argument
if (!is.na(args$dir_fastq)){
    dir_fastq <- normalizePath(args$dir_fastq)
} else {
    stop("Fastq directory must be provided. Use '--fastq <PATH>'.")
}
if (!is.na(args$dir_output)){
    dir_output <- normalizePath(args$dir_output)
} else {
    stop("Output directory must be provided. Use '--out <PATH>'.")
}
prefix <- args$prefix
trunc_len <- args$trunc_len
fastq_suffix_pattern <- args$fastq_suffix_pattern

## import rest of packages
library(tidyverse)
library(dada2)

## get fastq file sample IDs
if (args$gz){
    fastq_pattern <- paste0("^.+", fastq_suffix_pattern, "\\.fastq\\.gz$")
} else {
    fastq_pattern <- paste0("^.+", fastq_suffix_pattern, "\\.fastq$")
}
samples <- list.files(dir_fastq, pattern = fastq_pattern) %>%
    str_extract(paste0("^.+(?=", fastq_suffix_pattern, ")")) %>% unique %>% sort

## make fastq file names
suffix_r1 <- str_replace(fastq_suffix_pattern, "R\\\\d", "R1")
suffix_r2 <- str_replace(fastq_suffix_pattern, "R\\\\d", "R2")
fwd_reads <- paste0(dir_fastq, '/', samples, suffix_r1, ".fastq")
rvs_reads <- paste0(dir_fastq, '/', samples, suffix_r2, ".fastq")
if (args$gz){
    fwd_reads <- paste0(fwd_reads, ".gz")
    rvs_reads <- paste0(rvs_reads, ".gz")
}

## print(fastq_suffix_pattern)
## print(suffix_r1)
## print(fwd_reads)
## stop()

## generate file names for later
dada_dir <- paste0(dir_output, "/dada2/")
filt_dir <- paste0(dada_dir, '/', prefix)
filt_dir <- paste0(dada_dir, '/', prefix, "_trunc", trunc_len)
filtered_fwd_reads <- paste0(filt_dir, '/', samples, suffix_r1, ".filtered.trunc", trunc_len, ".fastq.gz")
filtered_rvs_reads <- paste0(filt_dir, '/', samples, suffix_r2, ".filtered.trunc", trunc_len, ".fastq.gz")

## ## visualisation to check read quality (Phred), not for ssh
## plotQualityProfile(fwd_reads)
## plotQualityProfile(rvs_reads)
filtered_out <- filterAndTrim(fwd_reads, filtered_fwd_reads, rvs_reads, filtered_rvs_reads,
                              maxEE = c(args$maxEE, args$maxEE), rm.phix = args$rm_phix,
                              minLen = trunc_len, truncLen = c(trunc_len, trunc_len))

## error model
err_fwd_reads <- learnErrors(filtered_fwd_reads, multithread = args$multithread)
err_rvs_reads <- learnErrors(filtered_rvs_reads, multithread = args$multithread)

## dereplication
derep_fwd <- derepFastq(filtered_fwd_reads, verbose = TRUE)
names(derep_fwd) <- samples
derep_rvs <- derepFastq(filtered_rvs_reads, verbose = TRUE)
names(derep_rvs) <- samples

## write derep to file
library(seqinr)
writeDerep <- function(dir_out, derep_fwd, derep_rvs = NULL){
    for (sample_id in names(derep_fwd)){
        print(sample_id)
        fwd <- derep_fwd[[sample_id]]
        rvs <- derep_rvs[[sample_id]]
        ## output files
        fa_fwd <- paste0(dir_out, "/", sample_id, ".derep.fwd.fasta")
        fa_rvs <- paste0(dir_out, "/", sample_id, ".derep.rvs.fasta")
        txt_map <- paste0(dir_out, "/", sample_id, ".derep.map.txt")
        txt_sum <- paste0(dir_out, "/", sample_id, ".derep.sum.txt")
        ## write seqs
        seqinr::write.fasta(sequences = fwd$unique %>% names %>% as.list,
                            names = paste0(1:length(fwd$unique), " n=", fwd$unique) %>% as.list,
                            as.string = TRUE, file.out = fa_fwd)
        seqinr::write.fasta(sequences = rvs$unique %>% names %>% as.list,
                            names = paste0(1:length(rvs$unique), " n=", rvs$unique) %>% as.list,
                            as.string = TRUE, file.out = fa_rvs)
        ## write mapping
        df <- data.frame(read = 1:length(fwd$map), fwd = fwd$map, rvs = rvs$map)
        write.table(df, file = txt_map, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
        ## write summary
        df.sum <- df %>%
            dplyr::group_by(fwd, rvs) %>%
            dplyr::summarise(reads = n()) %>%
            dplyr::arrange(desc(reads))
        write.table(df.sum, file = txt_sum, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
    }
}

dir_derep <- paste0(dir_output, "/derep/", prefix, "_trunc", trunc_len)
dir.create(dir_derep, showWarnings = FALSE, recursive = TRUE)
writeDerep(dir_derep, derep_fwd, derep_rvs)

