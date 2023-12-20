# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:46:07 2022

@author: Jenn Knapp
email: jknapp@uwaterloo.ca
"""
"""
Purpose: Organize and prep fastq files from an illumina sequencing run

Prep involves:
For all gzipped fastq files in a specified directory:
    Pair the reads, align them to a reference, index the bam.
Requires cutadapt, minimap2, samtools, and a reference genome (called sars-cov-2.fasta)
"""
#### Edit below for your system ####

data_path = "all_fastq/" # Directory with all fastq files to process
output_master = "Output/" # All outputs (sam, bam, bai) will end up here in separate directories for each sample
REF = "/home/jknapp/data/REF/sars-cov-2.fasta" # Path to reference wildtype sequence in fasta format

#### Edit above for your system ####

from os import listdir, mkdir, system
from os.path import isfile, join

import subprocess
import time
import os
import glob

if not os.path.isdir(data_path):
    os.mkdir(data_path)
    
if not os.path.isdir(output_master):
    os.mkdir(output_master)

samples = []
i = 1

# Check for fastq.gz files in data_path directory and add them to a list
all_fastq = glob.glob(data_path+"*.fastq.gz")
if all_fastq == []:
    raise FileNotFoundError("NO FASTQ.GZ FILES FOUND IN THE INPUT FOLDER!")
    
for read1 in all_fastq:
    if not "_R1" in os.path.basename(read1):
        # ignore R2
        continue
    else:
    # extract sample name
    # assume all input fastq follow the same pattern: LAB-####_*_R1_*.fastq.gz
    # each sample output is sorted into a separate folder with the name "LAB-####"
        base = os.path.basename(read1).split("_")[0]
        output_dir_name = base+'/'
    
    # If output directories do not exist, make them    
        output_dir = os.path.join(output_master, output_dir_name)
        if not os.path.isdir(output_dir):
            os.system("mkdir -p "+output_dir)
            continue
    
    # Find the corresponding R2
        read2 = data_path+os.path.basename(read1).replace('R1', 'R2')
        print("Aligning sample", base, ":", i, "out of", int(len(all_fastq) / 2)) 
        print("Using read1", read1, "and read2", read2)
    
    # Trim off adapters using cutadapt
        subprocess.run(['cutadapt', '-a', 'AGATCGGAAGAGC', '-A', 'AGATCGGAAGAGC', '-m', '1', '-o', read1.replace('.fastq.gz', 'prepped.fastq.gz'), '-p', read2.replace('.fastq.gz', 'prepped.fastq.gz'), read1, read2])	
        read1_prepped = read1.replace('.fastq.gz', 'prepped.fastq.gz')
        read2_prepped = read2.replace('.fastq.gz', 'prepped.fastq.gz')
  
    # Run minimap2 commands: minimap2 -ax sr 'sars-cov-2.fasta' ~/Data/path/all_fastq/sample1_1.fastq.gz ~/Data/path/all_fastq/sample1_2.fastq.gz > sample1.sam
      #   subprocess.run(['minimap2', '-ax', 'sr', 'sars-cov-2.fasta', 'read1_prepped', 'read2_prepped', '>', 'output_dir+base+.sam'])
        system('minimap2 -ax sr {0} {1} {2} > {3}{4}.sam'.format(REF, read1_prepped, read2_prepped, output_dir, base))   
    
    # Convert sam files to sorted bam files
        print("Converting sam file to sorted bam file for", base)
        subprocess.run(['samtools', 'sort', '-o', output_dir+base+'-sorted.bam', output_dir+base+'.sam'])
    
    # Index bam files    
        print("Indexing ", output_dir+base+'-sorted.bam')    
        subprocess.run(['samtools', 'index', output_dir+base+'-sorted.bam'])

    # Remove large intermediate files
        subprocess.run(['rm', read1_prepped])
        subprocess.run(['rm', read2_prepped])
        subprocess.run(['rm', output_dir+base+'.sam'])

    # Adding the path for the bam file and the sample name to a list    
        samples.append([output_dir+base+'-sorted.bam', base])
    
        i = i + 1
    
# Make samples.txt in parent directory - needed as input for ALCOV                
with open('samples_mm2.txt', 'w') as f:
    f.write('\n'.join(['\t'.join(line) for line in samples]))
   
subprocess.run(['sort','samples_mm2.txt', '-o', 'samples_mm2.txt']) 
