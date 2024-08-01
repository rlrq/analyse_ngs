#!/usr/bin/env python
# coding: utf-8

## adapted from Crispresso2_wrapper_v1.py

import argparse

parser = argparse.ArgumentParser(description="extract booking sheet from Excel file and convert to long tsv")
parser.add_argument("--amplicon_book", "--booking-excel", dest="excel",
                    required=True, help="path to booking Excel file")
parser.add_argument("--sheet", "--sheetname", "--name", dest="sheetname",
                    required=True, help="sheet name in Excel file")
parser.add_argument("--out", dest="fout",
                    required=True, help="path to output file")
args = parser.parse_args()

import pandas as pd

## parse file
amplicon_book = pd.read_excel(args.excel, sheet_name = args.sheetname)
amplicon_book_long = pd.wide_to_long(
    amplicon_book,
    stubnames = ['user', 'accession/gene', 'PCR_Product_fa', 'Guide'], i = ['fastq', 'Index1', 'Index2'],
    j = ' Set', suffix = r' \Set\w'
)

## remove all whitespaces 
df_obj  = amplicon_book_long.select_dtypes(['object'])
amplicon_book_long[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
amplicon_book_long = amplicon_book_long.reset_index()

## subset columns & remove null rows
amplicon_book_long = amplicon_book_long[['fastq', 'PCR_Product_fa', 'Guide']] 
amplicon_book_long = amplicon_book_long[amplicon_book_long.PCR_Product_fa.notnull()]
amplicon_book_long['PCR_Product_fa'] = amplicon_book_long['PCR_Product_fa'].str.replace('.fa$', '', regex=True)
amplicon_book_long['Guide'] = amplicon_book_long['Guide'].str.upper()

## write to file
amplicon_book_long.to_csv(args.fout, sep='\t')
