#!/usr/bin/python
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math
import subprocess
import pprint
import pandas as pd
import numpy as np



def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Convert BioHansel Kmer report to fasta sequence')
    parser.add_argument('--kmer_file', type=str, required=True, help='BioHansel Kmer output file')
    parser.add_argument('--scheme_map_file', type=str, required=True, help='BioHansel Scheme file of tiles to bases')
    parser.add_argument('--sample_name', type=str, required=False, help='Overwrite sample name in BioHansel results',default='')
    parser.add_argument('--out_file', type=str, required=False, help='If set, write results to file instead of screen',default='')
    return parser.parse_args()




def parse_scheme_lookup(scheme_lookup_file):
    scheme_lookup = pd.read_csv(scheme_lookup_file, sep='\t', index_col='tile_name',
                                names=['tile_name', 'position', 'base_called'], header=0)

    scheme = dict()
    for tile_name,row in scheme_lookup.iterrows():

        scheme[tile_name] = {
            'position': row['position'],
            'base_called': row['base_called'],
        }

    return scheme

def get_positions(scheme):
    positions = list()
    for tile in scheme:
        info = scheme[tile]
        position = int(info['position'])
        positions.append(position)

    positions = list(set(positions))
    positions.sort()
    return positions


def main():
    cmd_args = parse_args()

    in_kmer_file = cmd_args.kmer_file
    scheme_lookup_file = cmd_args.scheme_map_file
    sample_name = cmd_args.sample_name
    out_file = cmd_args.out_file



    kmer_info = pd.read_csv(in_kmer_file, sep='\t', index_col='kmername',
                            names=['kmername',
                                   'seq',
                                   'freq',
                                   'refposition',
                                   'subtype',
                                   'is_pos_kmer',
                                   'is_kmer_freq_okay',
                                   'file_path',
                                   'sample',
                                   'scheme',
                                   'scheme_version',
                                   'qc_status',
                                   'qc_message', ], header=0)

    scheme = parse_scheme_lookup(scheme_lookup_file)


    positions = get_positions(scheme)


    sequence = dict()
    for kmername, row in kmer_info.iterrows():
        if(len(sample_name) == 0):
            sample_name = row['sample']

        if kmername not in scheme:
            print("Error {} tile not in lookup table".format(kmername))
            continue
        else:
            position = scheme[kmername]['position']
            base = scheme[kmername]['base_called']

        if position not in sequence:
            sequence[position] = base
        else:
            if sequence[position] != base:
                sequence[position] = '-'

    seq = list()
    for position in positions:

        if position not in sequence:
            seq.append('-')
        else:
            seq.append(sequence[position])

    if len(out_file) > 0:
        fh = open(out_file,'w')
        fh.write(">{}\n{}".format(sample_name,"".join(seq)))
        fh.close()
    else:
        print(">{}\n{}".format(sample_name,"".join(seq)))


main()