#!#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(
    description="generate mutation frequency table from Crispresso output")
required_named_args = parser.add_argument_group("required named arguments")
required_named_args.add_argument("-c", "--crispresso", "--crispresso-dir", type=str, dest="dir_crispresso", required=True,
                                 help="Crispresso output directory")
required_named_args.add_argument("-e", "--excel-tsv", "-m", "--metadata", type=str, dest="metadata", required=True,
                                 help="tsv file of metadata")
required_named_args.add_argument("-o", "--out", type=str, dest="fout", required=True,
                                 help="output file")
required_named_args.add_argument("--fastq2sample", type=str, dest="f_fastq2sample", required=True,
                                 help = "path to fastq2sample.py (required)")
parser.add_argument("--sample-id-format", type=int, dest="sample_id_format", required=True, default=1,
                    help = "sample ID format; see fastq2sample.R for formats")
parser.add_argument("-p", "--gene-pattern", type=str, dest="gene_pattern", default = "[^_]+$",
                    help="gene regex; used with re.search to get gene ID from PCR_Product_fa column of metadata")
args = parser.parse_args()

dir_crispresso = args.dir_crispresso
gene_pattern = args.gene_pattern
f_fastq2sample = args.f_fastq2sample
sample_id_format = args.sample_id_format

import os
import re

from importlib.machinery import SourceFileLoader
f2s = SourceFileLoader("f2s", f_fastq2sample).load_module()
fastq2sample = f2s.fastq2sample(fmt = sample_id_format)

f_metadata = "/media/HDD3/rachelle/scripts/analyse_ngs/test/test4/metadata/20240724_Ploop.tsv"
f_fastq2sample = '/media/HDD3/rachelle/scripts/analyse_ngs/scripts/fastq2sample.py'
gene_pattern = "[^_]+$"
dir_crispresso = "/media/HDD3/rachelle/scripts/analyse_ngs/test/test4/crispresso/20240724_Ploop_trunc70"


with open(args.f_metadata, 'r') as f:
    header = f.readline().split('\t')
    i_fastq = header.index("fastq")
    i_pcr_product_fa = header.index("PCR_Product_fa")
    with open(args.fout, "w+") as fout:
        ## write header
        fout.write('\t'.join(["fastq", "pcr_product_fa", "gene", "site", "ref", "alt", "mut", "freq"]) + '\t')
        ## iterate through all samples
        for entry in f.readlines():
            entry_split = entry.split('\t')
            fastq = entry_split[i_fastq]
            pcr_product_fa = entry_split[i_pcr_product_fa]
            sample_id = fastq2sample(fastq)
            gene = re.search(gene_pattern, pcr_product_fa).group(0)
            f_pct_table = os.path.join(
                dir_crispresso, pcr_product_fa,
                f"CRISPResso_on_{sample_id}.{gene}.fwd_{sample_id}.{gene}.rvs",
                f"{sample_id}.{gene}_{pcr_product_fa}.{pcr_product_fa}.Quantification_window_nucleotide_percentage_table.txt")
            if not os.path.exists(f_pct_table):
                print( ("Could not find nucleotide percentage table for sample '{sample_id}'"
                        " with PCR product '{pcr_product_fa}'. Skipping.") )
                continue
            ## parse data
            with open(f_pct_table, 'r') as f:
                reference = f.readline().split('\t')[1:]
                data = {}
                for line in f:
                    split_line = line.split('\t')
                    data[split_line[0]] = [float(x) for x in split_line[1:]]
            ## convert from wide to long
            ## columns: Site, Ref, Alt, Mut, Freq
            for i, ref in enumerate(reference):
                for alt, freq in data.items():
                    fout.write('\t'.join(
                        map(str, [fastq, pcr_product_fa, gene, i+1, ref, alt, '-' if (ref==alt) else (ref+alt), freq[i]])) + '\n')

# ## write
# with open(fout, "w+") as f:
#     f.write('\t'.join(["Site", "Ref", "Alt", "Mut", "Freq"]) + '\t')
#     for entry in output:
#         f.write('\t'.join(map(str, entry)) + '\n')

exit 0

# fnames_pct_table = [f for f in os.listdir(dir_crispresso)
#                     if re.search("Quantification_window_nucleotide_percentage_table.txt$", f)]

# if len(fnames_pct_table) == 0:
#     print("XXX.Quantification_window_nucleotide_percentage_table.txt not found. Exiting.")
#     exit 1

# with open(os.path.join(dir_crispresso, fnames_pct_table[0]), 'r') as f:
#     reference = f.readline().split('\t')[1:]
#     data = {}
#     for line in f:
#         split_line = line.split('\t')
#         data[split_line[0]] = [float(x) for x in split_line[1:]]

# ## columns: Site, Ref, Alt, Mut, Freq
# output = []
# for i, ref in enumerate(reference):
#     for alt, freq in data.items():
#         output.append([i+1, ref, alt, '-' if (ref==alt) else (ref+alt), freq[i]])

# ## write
# with open(fout, "w+") as f:
#     f.write('\t'.join(["Site", "Ref", "Alt", "Mut", "Freq"]) + '\t')
#     for entry in output:
#         f.write('\t'.join(entry) + '\n')

# exit 0
