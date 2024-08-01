#!/usr/bin/env R

fastq2sample <- function(format = 1){
    if (format==1){
        output <- function(fastq_id){
            return(paste0(fastq_id, '_', fastq_id) %>% stringr::str_replace("_S0", "_S"))
        }
    }
    return(output)
}
