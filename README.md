# CoV Breakdown

Analysis of SARS-CoV-2 mutations. The tool can be used for:

* Converting nucleotide and amino acid mutations for SARS-CoV-2 such as those found on https://covariants.org/variants/S.N501
* Determining the frequency of mutations of interest in BAM files. Currently looks for N501Y and A28271- from https://www.medrxiv.org/content/10.1101/2021.02.22.21252041v1.full.pdf
* Plotting the depth for each ARTIC amplicon (https://github.com/artic-network/artic-ncov2019/tree/master/primer\_schemes/nCoV-2019/V3)
* Comparing amplicon GC content with its read depth (as a measure of degredation)

## Installing

Clone the repository and run

`pip install .`

This will install the Python library and the CLI.

Note: `find_mutants` also requires pysam and seaborn

## Usage example

Converting mutation names:

```
$ cov_breakdown nt A23063T
A23063T causes S:N501Y
$ cov_breakdown aa S:E484K
G23012A causes S:E484K
```

### Finding mutations in BAM file:

```
cov_breakdown find_mutants reads.bam
```

Finding mutations in BAM files for multiple samples:

```
cov_breakdown find_mutants samples.txt
```

Where `samples.txt` looks like:

```
reads1.bam	Sample 1 name
reads2.bam	Sample 2 name
...
```

Running `find_mutants` will print the number of reads with and without each mutation in each sample and then generate a heatmap showing the frequencies for all samples.

You can also specify a custom mutations file:

```
cov_breakdown find_mutants samples.txt mutations.txt
```

Where `mutations.txt` looks like:

```
S:N501Y
G23012A
...
```

### Getting the read depth for each amplicon

```
cov_breakdown amplicon_coverage reads.bam
```

or

```
cov_breakdown amplicon_coverage samples.txt
```

### Plotting amplicon GC content against amplicon depth

```
cov_breakdown gc_depth reads.bam
```

or

```
cov_breakdown gc_depth samples.txt
```
