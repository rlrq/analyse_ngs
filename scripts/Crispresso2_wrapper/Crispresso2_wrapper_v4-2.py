#!/usr/bin/env python
# coding: utf-8

import argparse
import subprocess
import glob
import os
import pandas as pd
from multiprocessing import Pool
from Bio import SeqIO
import traceback
import sys
import logging

# Create the parser and add arguments
parser = argparse.ArgumentParser(description='Run CRISPResso analysis.')
parser.add_argument('--fastq_dir', required=True, help='Directory containing FASTQ files')
parser.add_argument('--amplicon_fasta_dir', required=True, help='Directory containing amplicon FASTA files')
parser.add_argument('--output_dir', default='./Crispresso_output', help='Output directory for CRISPResso (default: ./Crispresso_output)')
parser.add_argument('--verbose', type=str, choices=['T', 'F'], default='T', help='Verbose output (T/F, default: T)')
## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
## (removed required=True for amplicon_book & sheet)
parser.add_argument('--amplicon_book', help='Amplicon booking list in xlsx format')
parser.add_argument('--sheet', help='Sheet name in the Amplicon booking list')
parser.add_argument('--excel-tsv', '--amplicon_tsv', help='Amplicon booking sheet in long TSV format', dest = "amplicon_tsv")
parser.add_argument('--fastq_name_format', type=int, help='Fastq file name format (default=2)', default=2)
# parser.add_argument('--editor', type=str, help='Editor. Valid values: Cas9, ABE (default=Cas9)',
#                     default="Cas9", dest="editor", choices=["Cas9", "ABE"])
parser.add_argument('--PCR_Product_fa_gene_pattern', type=str, default='([^_]+)$', dest="gene_pattern",
                    help="pattern to extract gene from PCR_Product_fa column (default='([^_]+)$')")

# Parse the arguments
args, args_free = parser.parse_known_args()

# Use the arguments in the script
fastq_dir = args.fastq_dir
amplicon_fasta_dir = args.amplicon_fasta_dir
output_dir = args.output_dir
amplicon_book_path = args.amplicon_book
sheet_name = args.sheet
verbose = args.verbose == 'T'
fastq_name_format = args.fastq_name_format
gene_pattern = args.gene_pattern
# editor = args.editor

## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
amplicon_tsv_path = args.amplicon_tsv
if (all(map(lambda x:x is None, [amplicon_tsv_path, amplicon_book_path, sheet_name]))
    or amplicon_tsv_path is not None and not (amplicon_book_path is None and sheet_name is None)):
    print("Either '--amplicon_tsv <path>' OR '--amplicon_book <path> --sheet <sheet name>' is required.")
    sys.exit()
elif amplicon_tsv_path is None and (amplicon_book_path is None or sheet_name is None):
    print("--amplicon_book must be used with --sheet.")

## make output directory
if not os.path.exists(output_dir):
    os.mkdir(output_dir) 
else:
    ## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
    print("Output directory already exists, files with same name will be overwritten.")
    # print('Output directory already exists, use --output_dir to change the output directory')
    # sys.exit()

# Set up logging
log_file = os.path.join(output_dir, 'crispresso_wrapper.log')
logging.basicConfig(level=logging.INFO if verbose else logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.FileHandler(log_file),
                        logging.StreamHandler()
                    ])

# Function to log messages both to file and console
def log_message(message, level=logging.INFO):
    if verbose:
        logging.log(level, message)

log_message(f"Starting CRISPResso analysis with arguments: {args}")


## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
## parse amplicon book if amplicon_book_path is set
if amplicon_book_path and sheet_name:
    amplicon_book = pd.read_excel(amplicon_book_path, sheet_name=sheet_name)
    amplicon_book_long = pd.wide_to_long(
        amplicon_book, stubnames=['user', 'accession/gene', 'PCR_Product_fa', 'Guide'], 
        i=['fastq', 'Index1', 'Index2'],
        j=' Set', suffix=r' \Set\w')
else:
    amplicon_book_long = pd.read_csv(amplicon_tsv_path, sep='\t')

# Remove all whitespaces 
df_obj = amplicon_book_long.select_dtypes(['object'])
amplicon_book_long[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
amplicon_book_long = amplicon_book_long.reset_index()

## reformat cells
amplicon_book_long = amplicon_book_long[['fastq', 'PCR_Product_fa', 'Guide', 'accession/gene']] 
amplicon_book_long = amplicon_book_long[amplicon_book_long.PCR_Product_fa.notnull()]
amplicon_book_long['PCR_Product_fa'] = amplicon_book_long['PCR_Product_fa'].str.replace('.fa$', '', regex=True)
amplicon_book_long['Guide'] = amplicon_book_long['Guide'].astype(str)
amplicon_book_long['Guide'] = amplicon_book_long['Guide'].str.upper()

amplicon_book_long['Guide'] = amplicon_book_long['Guide'].str.upper().str.strip()

## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
## make column for gene
amplicon_book_long["Gene_ID"] = amplicon_book_long["PCR_Product_fa"].str.extract(gene_pattern)

actual_files = set(os.path.basename(f) for f in glob.glob(amplicon_fasta_dir + '/*.fa'))
valid_files_mask = amplicon_book_long['PCR_Product_fa'].apply(lambda x: x + '.fa' in actual_files)
invalid_rows = amplicon_book_long[~valid_files_mask]
amplicon_book_long = amplicon_book_long[valid_files_mask]

if not invalid_rows.empty:
    log_message(f'These fasta files are missing: {set(invalid_rows["PCR_Product_fa"])}', logging.WARNING)

def get_amplicon_seq(file_path):
    try:
        seq_record = SeqIO.read(file_path, "fasta")
        return str(seq_record.seq)
    except ValueError as e:
        raise ValueError(f"Error reading FASTA file {file_path}: {str(e)}")

def get_prefix_fastq(s_value):
    number = int(s_value.strip("S"))
    return f"{s_value}_{s_value}" if number >= 10 else f"{s_value}_S{number}"

## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
def make_fastq_filename(s_value, gene, fmt = 1):
    if fmt == 1:
        prefix = f"{s_value}_S{int(s_value.strip('S'))}"
        return prefix, prefix + "_L001_R1_001.fastq.gz", prefix + "_L001_R2_001.fastq.gz"
    elif fmt == 2:
        prefix = f"{s_value}_S{int(s_value.strip('S'))}.{gene}"
        return prefix, prefix + ".fwd.fastq.gz", prefix + ".rvs.fastq.gz"

def run_crispresso(row):
    ## ANALYSE_NGS MOD: modifications for analyse_ngs output compatibility
    fq_file, amplicon_fasta, guide, gene_id = row['fastq'], row['PCR_Product_fa'], row['Guide'], row['Gene_ID']
    amplicon_file = os.path.join(amplicon_fasta_dir, amplicon_fasta + '.fa')
    fastq_prefix, fastq_fwd, fastq_rvs = make_fastq_filename(fq_file, gene_id, fmt = fastq_name_format)
    # fastq_prefix = get_prefix_fastq(fq_file)
    log_file = os.path.join(output_dir, amplicon_fasta, fastq_prefix + '_' + amplicon_fasta + '_error_log.txt')
    
    log_message(f"Processing sample: {amplicon_fasta}")
    
    try:
        amplicon_seq = get_amplicon_seq(amplicon_file)
        
        crispresso_command = [
            'CRISPResso',
            '--fastq_r1', os.path.join(fastq_dir, fastq_fwd),
            '--fastq_r2', os.path.join(fastq_dir, fastq_rvs),
            '--amplicon_seq', amplicon_seq,
            '--keep_intermediate',
            '--amplicon_name', amplicon_fasta,
            '--bam_output',
            '-o', os.path.join(output_dir, amplicon_fasta),
            '--file_prefix', fastq_prefix + '_' + amplicon_fasta,
        ]
        
        if not pd.isna(guide):
            crispresso_command.extend(['--guide_seq', guide])

        ## add arbitrary arguments
        crispresso_command.extend(args_free)
        
        log_message(f"Running CRISPResso command for {amplicon_fasta}")
        crispresso_process = subprocess.run(crispresso_command, capture_output=True, text=True)
        if crispresso_process.returncode == 0:
            log_message(f"CRISPResso command was successful for {amplicon_fasta}")
            return f"CRISPResso command was successful for {amplicon_fasta}"
        else:
            error_message = f"CRISPResso command failed for {amplicon_fasta}! Error: {crispresso_process.stderr}"
            log_message(error_message, logging.ERROR)
            with open(log_file, 'a') as f:
                f.write(error_message + "\n")
            return error_message
    except Exception as e:
        error_message = f"An error occurred for {amplicon_fasta}: {str(e)}"
        log_message(error_message, logging.ERROR)
        with open(log_file, 'a') as f:
            f.write(error_message + "\n" + traceback.format_exc() + "\n")
        return error_message

log_message("Starting to process samples")
rows_as_dicts = amplicon_book_long.to_dict('records')
pool = Pool(processes=40) 
results = pool.map(run_crispresso, rows_as_dicts)

# Close the pool and wait for all processes to finish
pool.close()
pool.join()

log_message("Finished processing all samples")

# Print results
for result in results:
    log_message(result)

# The pattern to match files deep in the directory structure.
editing_freq_pattern = '**/*quantification_of_editing_frequency.txt'

# Full pattern path combining base directory and the file pattern.
full_editing_freq_pattern = os.path.join(output_dir, editing_freq_pattern)

# Perform the recursive search to get a list of all matching file paths.
editing_file_paths = glob.glob(full_editing_freq_pattern, recursive=True)

# Create a container for the DataFrames
dfs = []
for file_path in editing_file_paths:
    try:
        # Read the current file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t')
        
        # Add an identifier for each file
        current_df['Source_File'] = os.path.basename(file_path)
        
        # Append the current DataFrame to the container
        dfs.append(current_df)
    except Exception as e:
        log_message(f"Could not read file {file_path}. Error: {e}", logging.ERROR)
        
if dfs:
    final_df = pd.concat(dfs, ignore_index=True) 
    # Regular expression pattern to match strings like "S01", "S02", ..., "S95", etc.
    pattern = r'(?<!\d)(S\d{2})'

    # Extract the pattern from the 'Source_File' column and create a new column 'Fastq' with the result
    final_df['Fastq'] = final_df['Source_File'].str.extract(pattern, expand=False)
    final_df['Insertions%'] = 100 * final_df['Insertions'] / final_df['Reads_aligned']
    final_df['Deletions%'] = 100 * final_df['Deletions'] / final_df['Reads_aligned']
    final_df['Substitutions%'] = 100 * final_df['Substitutions'] / final_df['Reads_aligned']
    
    # Merge the accession/gene info with amplicon book
    amplicon_book_simple = amplicon_book_long.rename(columns={"fastq": 'Fastq', "PCR_Product_fa": "Amplicon"})
    amplicon_book_simple = amplicon_book_simple[['Fastq', 'Amplicon', 'accession/gene']]
    # Performing the left join
    final_df = pd.merge(final_df, amplicon_book_simple, on=["Fastq", "Amplicon"], how='left')

    # Define the desired column order
    desired_cols_order = [
        'Fastq', 'accession/gene', 'Amplicon', 'Unmodified%', 'Modified%', 
        'Insertions%', 'Deletions%', 'Substitutions%'
    ]

    # Get the remaining columns not specified in the desired order
    remaining_cols = [col for col in final_df.columns if col not in desired_cols_order]

    # Combine the two lists to get the full rearranged order
    cols_order = desired_cols_order + remaining_cols

    # Apply the new column order to the DataFrame
    final_df = final_df[cols_order]
    final_df.to_csv(os.path.join(output_dir, 'summary_editing_frequency.csv'), index=False)
    log_message("Summary of editing frequency saved to 'summary_editing_frequency.csv'")
else:
    log_message("No result files were read successfully.", logging.WARNING)

error_pattern = '**/*error_log.txt'

# Full pattern path combining base directory and the file pattern.
full_error_pattern = os.path.join(output_dir, error_pattern)

# Perform the recursive search to get a list of all matching file paths.
error_file_paths = glob.glob(full_error_pattern, recursive=True)

# Initialize an empty list to store the row data for our future DataFrame.
data = []

# Function to find the line starting with 'ERROR'.
def find_error_line(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ERROR'):
                return line.strip()  # Return the line without leading/trailing white spaces.
    return None  # Return None if no such line exists.

# Iterate over all the file paths.
for file_path in error_file_paths:
    try:
        # Call the function to get the error line.
        error_line = find_error_line(file_path)

        # Proceed only if an 'ERROR' line was found.
        if error_line:
            # Extract the filename from the full file path.
            base_name = os.path.basename(file_path)
            # Remove the specific part of the file name.
            error_name = base_name.replace('_error_log.txt', '')

            # Append the data to our rows list.
            data.append([error_name, error_line])
    except Exception as e:
        log_message(f"Failed to process file {file_path} due to {e}", logging.ERROR)

# Create a DataFrame from the accumulated data.
if data == []:
    log_message('No errors found.')
else:
    error_df = pd.DataFrame(data, columns=['File_Name', 'Error_Line'])
    error_df.to_csv(os.path.join(output_dir, 'summary_error_log.csv'), index=False)
    log_message("Summary of errors saved to 'summary_error_log.csv'")

log_message("CRISPResso analysis completed")
