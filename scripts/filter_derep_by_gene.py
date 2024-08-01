#!/usr/bin/python3

import sys
sys.path.append("/mnt/chaelab/rachelle/src")

from Bio import SeqIO
import data_manip

import argparse

parser = argparse.ArgumentParser(
    description="generate 'pattern' file for seqkit for demultiplex step of analyse_ngs")
required_named_args = parser.add_argument_group("required named arguments")
required_named_args.add_argument("-f", "--fasta", type=str, dest="fasta", required=True,
                                 help="fasta file of sequences to filter")
required_named_args.add_argument("-m", "--map", type=str, dest="map", required=True,
                                 help="tsv file of sequence-identity map")
required_named_args.add_argument("-g", "--gene", type=str, dest="gene", required=True,
                                 help="gene to filter reads for")
required_named_args.add_argument("-o", "--out", type=str, dest="fout", required=True,
                                 help="output file")
args = parser.parse_args()

def get_homologue_seq_id(f_seq_id, gene):
    get, data = data_manip.parse_get_data(f_seq_id, parse_num = False)
    output = [get(x, "qseqid") for x in data if (gene in get(x, "gene").split(','))]
    return set(output)

def filter_seqs(fasta, seqids):
    seqids = set(seqids)
    output = []
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in seqids:
            output.append(record)
    return output

def write_seqs_txt(fout, seq_records):
    with open(fout, "w+") as f:
        for record in seq_records:
            f.write(str(record.seq) + '\n')
    return

seqids_to_keep = get_homologue_seq_id(args.map, args.gene)
seqs_to_keep = filter_seqs(args.fasta, seqids_to_keep)
write_seqs_txt(args.fout, seqs_to_keep)
