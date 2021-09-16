# Alcov

Abundance learning for SARS-CoV-2 variants. The primary purpose of the tool is:

* Estimating abundace of variants of concern from wastewater sequencing data

You can read more about how Alcov works in the preprint, __[Alcov: Estimating Variant of Concern Abundance from SARS-CoV-2 Wastewater Sequencing Data](https://www.medrxiv.org/content/10.1101/2021.06.03.21258306v1)__

The tool can also be used for:

* Converting nucleotide and amino acid mutations for SARS-CoV-2 such as those found on https://covariants.org/variants/S.N501
* Determining the frequency of mutations of interest in BAM files
* Plotting the depth for each ARTIC amplicon (https://github.com/artic-network/artic-ncov2019/tree/master/primer\_schemes/nCoV-2019/V3)
* Comparing amplicon GC content with its read depth (as a measure of degredation)

The tool is under active development. If you have questions or issues, please open an issue on GitHub or email me (email in setup.py).

## Installing

The latest release can be downloaded from PyPI

`pip install alcov`

This will install the Python library and the CLI.

To install the development version, clone the repository and run

`pip install .`

## Usage example

### Estimating relative abundance of variants of concern:

```
alcov find_lineages reads.bam
```

Optionally look for non-voc lineages and change minimum read depth (default 40)

```
alcov find_lineages --only_vocs=False --min_depth=5 reads.bam
```

### Converting mutation names:

```
$ alcov nt A23063T
A23063T causes S:N501Y
$ alcov aa S:E484K
G23012A causes S:E484K
```

### Finding mutations in BAM file:

```
alcov find_mutants reads.bam
```

Finding mutations in BAM files for multiple samples:

```
alcov find_mutants samples.txt
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
alcov find_mutants samples.txt mutations.txt
```

Where `mutations.txt` looks like:

```
S:N501Y
G23012A
...
```

### Getting the read depth for each amplicon

```
alcov amplicon_coverage reads.bam
```

or

```
alcov amplicon_coverage samples.txt
```

### Plotting amplicon GC content against amplicon depth

```
alcov gc_depth reads.bam
```

or

```
alcov gc_depth samples.txt
```
