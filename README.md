# tRNA Modification Analysis and Predictions (tMAP)

getMismatches.py is an additional tool to extend the usage of the standard tRAX analysis pipeline. It can be used to determine potential modifcaiton-induced misincorporations from small RNA sequencing data and in some cases provide a probability of associated with a specfic modification.

## About

[tRAX](https://github.com/UCSC-LoweLab/tRAX) is a tool often used for analyzing tRNA-seq data. While it generates a comprehensive set of results, it does not provide a way to specifically assess modification-induced misincorporations present in the data, or provide these analyses on non-tRNAs. Here we provide a tool to expand the usage of [tRAX](https://github.com/UCSC-LoweLab/tRAX) to generate a probability of a base harboring an RNA modification, (optionally) predict modifications based on RT-based signatures, and (optionally) provide nucleotide-resolution coverage data for user-defined regions.

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f getmismatch.yml
```

## Input files

getmismatch will work with any tRAX generated coverage file

### Metadata

You will need to provide the meta-data associated with the samples if you want to group samples based on specific experimental conditions. To do this you can provide a [tRAX samples file](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) whitespace/.csv/.tsv file (`-s/--samples`) containing the sample names, sample groups, and path to fastq files. An example samples file is shown below:

```tsv
sample1 a /path/to/sample1.fastq.gz
sample2 a /path/to/sample2.fastq.gz
sample3 b /path/to/sample3.fastq.gz
sample4 b /path/to/sample4.fastq.gz
```

## Usage

### Activating the environment

```bash
conda activate getmismatch
```

### Determining potential sites of base modifications present in tRNAs

The getMismatches.py script can be used to identify potential sites of base modification present in tRNA sequencing data.

```bash
python getMismatches.py --cov <tRAX_coverage_file> --o <output_file_prefix> [(optional) --alpha <mismatch_pseudocounts> --beta <reference_pseudocounts> --editfrac <minumum_edit_fraction> --minreads <minimum_read_coverage> --multimapcutoff <percent_multimap_coverage> --exppairs <pairs.txt> --samples <samples.txt> --org <organism> --predict]
```
* `--cov` is the path to the tRAX coverage file you want to analyse modifications from; runname-coverage.txt.
* `--o`  is the prefix for the output files.
* (optional flags)
* `--alpha` pseudocounts to add to mismatch counts for beta function; default=0.
* `--beta` pseudocounts to add to reference counts for beta function; default=0.
* `--editfrac` minimum fraction of mismatched reads for beta function; default=0.05.
* `--minreads` is the minumum number of read coverage of a base required to go into analysis; default=20.
* `--positions` tRNA positions of interest. Can list multiple positions e.g. 9 26 32 34...
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--exppairs` pairs.txt file from tRAX to perform pairwise comaprisons.
* `--samples` samples.txt file from tRAX, or custom made as previously described to group biological replicates.
* `--multimapcutoff` Percent of multimapping reads covering a site to mark as artefact; default=20.
* `--predict` (in development) predict which base modification is present based on misincorporation signature. Model built using tRNA sequencing (OTTR-seq) data from various eukaryotic model organisms.
