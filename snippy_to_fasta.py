#!/usr/bin/python
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math
from Bio import SeqIO
import pandas as pd
import numpy as np

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Use SNIPPY output files with a list of target sites to produce a fasta sequence.'
                                        ' All SNIPPY files need to be from the same SNIPPY run')
    parser.add_argument('--align_fa', type=str, required=True, help='SNIPPY *.align.fasta')
    parser.add_argument('--tab', type=str, required=True, help='SNIPPY Tab output')
    parser.add_argument('--posiitons', type=str, required=True, help='Two column table of CHROM,POS')
    parser.add_argument('--sample_name', type=str, required=False, help='Overwrite sample name in BioHansel results',
                        default='')
    parser.add_argument('--out_file', type=str, required=False, help='If set, write results to file instead of screen',default='')
    return parser.parse_args()



def read_snippy_tab(snippy_tab_file):
    snippy = pd.read_csv(snippy_tab_file, sep='\t', index_col='POS',
                                names=[
                                    'CHROM',
                                    'POS',
                                    'TYPE',
                                    'REF',
                                    'ALT',
                                    'EVIDENCE',
                                    'FTYPE',
                                    'STRAND',
                                    'NT_POS',
                                    'AA_POS',
                                    'EFFECT',
                                    'LOCUS_TAG',
                                    'GENE',
                                    'PRODUCT',
                                ], header=0)
    snps = dict()
    for pos,row in snippy.iterrows():
        if not row['CHROM'] in snps:
            snps[row['CHROM']] = dict()
        if row['TYPE'] != 'snp':
            continue
        snps[row['CHROM']][pos] = row['ALT']

    return snps

def read_positions_file(positions_file):
    pos_file = pd.read_csv(positions_file, sep='\t', header=0)
    positions = dict()
    for i,row in pos_file.iterrows():

        if not row['CHROM'] in positions:
            positions[row['CHROM']] = list()

        pos = row['POS']
        positions[row['CHROM']].append(pos)

    for chr in positions:
        positions[chr] = list(set(positions[chr]))
        positions[chr].sort()
    return positions



def main():
    cmd_args = parse_args()
    in_fasta = cmd_args.align_fa
    snippy_tab_file = cmd_args.tab
    positions_file = cmd_args.posiitons
    out_file = cmd_args.out_file
    sample_name = cmd_args.sample_name

    snps = read_snippy_tab(snippy_tab_file)
    positions = read_positions_file(positions_file)



    for seq_record in SeqIO.parse(in_fasta, "fasta"):
        if(len(sample_name) == 0):
            sample_name = seq_record.id

        if seq_record.id not in snps:
            print("Error your fasta file has sequence id: {} which is not in your snippy tab file".format(seq_record.id))
            continue
        if seq_record.id not in positions:
            print("Error your fasta file has sequence id: {} which is not in your positions file".format(seq_record.id))
            continue

        seq_len = len(seq_record.seq)
        seq = list()

        for pos in positions[seq_record.id]:
            if pos > seq_len or pos < 1:
                print("Error you have specified the position: {} for {} record which out of bounds".format(pos,seq_record.id))
            if pos in snps[seq_record.id]:
                base = snps[seq_record.id][pos]
            else:
                base = seq_record.seq[pos-1]
            seq.append(base)

        if len(out_file) > 0:
            fh = open(out_file, 'w')
            fh.write(">{}\n{}".format(sample_name, "".join(seq)))
            fh.close()
        else:
            print(">{}\n{}".format(sample_name,"".join(seq)))

main()