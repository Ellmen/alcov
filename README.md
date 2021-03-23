# CoV Breakdown

Analysis of SARS-CoV-2 mutations. The tool can be used for:

* Converting nucleotide and amino acid mutations for SARS-CoV-2 such as those found on https://covariants.org/variants/S.N501
* Determining the frequency of mutations of interest in BAM files. Currently looks for N501Y and A28271- from https://www.medrxiv.org/content/10.1101/2021.02.22.21252041v1.full.pdf

## Installing

Clone the repository and run

`pip install .`

This will install the Python library and the CLI.

## Usage example

Converting mutation names:

```
$ cov_breakdown nt A23063T
A23063T causes S:N501Y
$ cov_breakdown aa S:E484K
G23012A causes S:E484K
```

Finding mutations in BAM file:

```
cov_breakdown find_mutants reads.bam
```
