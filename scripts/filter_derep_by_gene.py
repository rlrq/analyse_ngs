#!/usr/bin/env python3

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

from Bio import SeqIO

def make_custom_get(header, parse_num = True):
    def get_col(colname, data):
        return [get_col_in_row(x, colname) for x in data]
    def get_col_in_row(row, colname):
        if not colname in header:
            print(f"Column '{colname}' is not found in headers {header}")
        output = row[header.index(colname)]
        if isinstance(output, (list, tuple)):
            if [str(x).isdigit() for x in output].count(True) == len(output):
                output = [int(x) for x in output]
            elif [str(x).replace('.','',1).replace('-','',1).isdigit() \
                  for x in output].count(True) == len(output):
                output = [float(x) for x in output]
            return output
        else:
            return output if (parse_num is False
                              or (not str(output).replace('.','',1).replace('-','',1).isdigit())) else \
                float(output) if not str(output).isdigit() else int(output)
    def helper(data = None, *colnames, get_cols = False, suppress_print = False, ncol = False,
               return_list = False):
        if get_cols:
            return header
        if ncol:
            return len(header)
        if len(data) == 0:
            if not suppress_print:
                print("No data found; returning empty list")
            return []
        if isinstance(data[0], (list, tuple)):
            output = [get_col(colname, data) for colname in colnames]
            return output[0] if len(output) == 1 else [[output[r][c] for r in range(len(output))] \
                                                       for c in range(len(output[0]))]
        else:
            output = [get_col_in_row(data, colname) for colname in colnames]
            return output[0] if (len(output) == 1 and not return_list) else output
    return helper

def parse_get_data(fname, delim = '\t', detect = True, parse_num = True):
    if detect:
        ext = fname.split('.')[-1]
        if ext == "csv":
            delim = ','
        elif ext == "tsv":
            delim = '\t'
    # data = [(''.join([c for c in line if c != '\n'])).split(delim)
    #         for line in open(fname, 'r').readlines()]
    data = [line.split(delim) for line in splitlines(fname)]
    get = make_custom_get(data[0], parse_num = parse_num)
    return get, data[1:]

def splitlines(fname, ignore_empty_lines = True):
    data = open(fname, 'r').read().split('\n')
    if ignore_empty_lines:
        return [line for line in data if line]
    else:
        return data

def get_homologue_seq_id(f_seq_id, gene):
    get, data = parse_get_data(f_seq_id, parse_num = False)
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
