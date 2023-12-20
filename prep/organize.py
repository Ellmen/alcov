from os import listdir, mkdir, system
from os.path import isfile, join

from shutil import copy

"""
For all gzipped fastq files in a specified directory:
Pair the reads, align them to a reference, index the bam.
Requires SeqPrep, bwa, samtools, and a reference genome (called sars-cov-2.fasta)
"""

data_path = './data/raw' # Reads directory (specify)
samples_path = './data/samples' # Processed (specify)

fns = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.fastq.gz')]
fns.sort()
samples = []

# Relies on pairs being alphabetically adjacent
pairs = []
for i in range(int(len(fns) / 2)):
    sample_path = join(samples_path, 'sample{}'.format(i+1))
    try:
        mkdir(sample_path)
    except:
        print('{} already exists'.format(sample_path))
    read1 = fns[i*2]
    read2 = fns[i*2+1]
    name = read1[:read1.find('_')]
    copy(join(data_path, read1), join(sample_path, read1))
    copy(join(data_path, read2), join(sample_path, read2))
    system('ls {}'.format(sample_path))
    print('Prepping reads')
    cmd = 'SeqPrep -f {0}/{1} -r {0}/{2} -1 {0}/prepped1.fastq.gz -2 {0}/prepped2.fastq.gz -s {0}/merged.fastq.gz'.format(sample_path,read1,read2)
    system(cmd)
    print('Aligning reads')
    merged = '{}/merged'.format(sample_path)
    system('bwa aln sars-cov-2.fasta {0}.fastq.gz > {0}.sai'.format(merged))
    system('bwa samse sars-cov-2.fasta {0}.sai {0}.fastq.gz > {0}.sam'.format(merged))
    system('samtools view -b {0}.sam > {0}.bam'.format(merged))
    system('samtools sort {0}.bam > {0}.sorted.bam'.format(merged))
    system('samtools index {0}.sorted.bam'.format(merged))
    samples.append(['sample{}/merged.sorted.bam'.format(i+1), name])

with open('{}/samples.txt'.format(samples_path), 'w') as f:
    f.write('\n'.join(['\t'.join(line) for line in samples]))
