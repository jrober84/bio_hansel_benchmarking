#!/usr/bin/python
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math
import subprocess
import pprint
import pandas as pd
import numpy as np



def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Tool for creating bash scripts for testing SNIPPY and biohansel accross difference coverage thresholds')
    parser.add_argument('--job_dir', type=str, required=True, help='Directory to write job files to')
    parser.add_argument('--sample_data', type=str, required=True, help='Sample information')
    parser.add_argument('--threads', type=int, required=False, help='Num threads for tools to use',default=1)
    return parser.parse_args()



def main():
    cmd_argments = parse_args()
    sample_file = cmd_argments.sample_data
    job_outdir = cmd_argments.job_dir
    threads = cmd_argments.threads
    sample_info = pd.read_csv(sample_file, sep='\t', index_col='sample_id',
                              names=['sample_id','outdir','snippy_reference','bio_hansel_scheme_file'],header=0)



    for sample_id, row in sample_info.iterrows():
        outdir = row['outdir']
        fastq_out_dir = os.path.join(outdir,"fastq")
        snippy_out = os.path.join(outdir,"snippy")
        bio_out = os.path.join(outdir,"bio_hansel")
        bio_hansel_scheme = row['bio_hansel_scheme_file']

        reference_file = row['snippy_reference']
        R1= os.path.join(fastq_out_dir,"{}_1.fastq".format(sample_id))
        R2 = os.path.join(fastq_out_dir,"{}_2.fastq".format(sample_id))
        snippy_coverages = [3,6,8]
        bio_coverages = [3, 6, 8]
        bash_string = "#!/bin/sh\n"

        #Get Data from NCBI
        bash_string = bash_string + "fastq-dump --split-files --outdir {} {} &&\n".format(fastq_out_dir, sample_id)

        for mincov in snippy_coverages:
            snippy_dir = os.path.join(snippy_out,"snippy_mincov_{}_{}".format(mincov,sample_id))
            bash_string+= "snippy --prefix {} --cpus {} --outdir {} --ref {} --R1 {} --R2 {} --cleanup --quiet --mincov {} &&\n".format(sample_id,threads,snippy_dir,reference_file,R1,R2,mincov)

        for mincov in bio_coverages:
            bio_dir = os.path.join(bio_out, "biohansel_mincov_{}_{}".format(mincov,sample_id))
            bash_string += "mkdir {} &&\n".format(bio_dir)
            summary_file = os.path.join(bio_dir, "{}.summary.txt".format(sample_id))
            simple_file = os.path.join(bio_dir, "{}.simple.txt".format(sample_id))
            kmer_file = os.path.join(bio_dir, "{}.kmer.txt".format(sample_id))
            bash_string += "hansel --min-kmer-freq {} -t {} -p {} {} -s {} -o {} -O {} -S {} &&\n".format(
                mincov,threads,R1,R2,bio_hansel_scheme,summary_file,kmer_file,simple_file)

        bash_string += "rm {} &&\n".format(R1)
        bash_string += "rm {} ;\n".format(R2)
        target = open(os.path.join(job_outdir,"{}.sh".format(sample_id)), 'w')
        target.write(bash_string)
        target.close()






# call main function
if __name__ == '__main__':
	main()
