#!/usr/bin/env python3

def fastq2sample(fmt):
    if fmt == 1:
        output = lambda fastq_id:(fastq_id + '_' + fastq_id.replace("S0", "S"))
    return output
