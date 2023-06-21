# tRNA Modification Analysis and Predictions (tMAP)

tMAP is an additional tool to extend the usage of the standard tRAX analysis pipeline. It can be used to determine potential modifcaiton-induced misincorporations from small RNA sequencing data and in some cases provide a probability of associated with a specfic modification.

## About

[tRAX](https://github.com/UCSC-LoweLab/tRAX) is a tool often used for analyzing tRNA-seq data. While it generates a comprehensive set of results, it does not provide a way to specifically assess modification-induced misincorporations present in the data, or provide these analyses on non-tRNAs. Here we provide a tool to expand the usage of [tRAX](https://github.com/UCSC-LoweLab/tRAX) to generate a probability of a base harboring an RNA modification, (optionally) predict modifications based on RT-based signatures, and (optionally) provide nucleotide-resolution coverage data for user-defined regions.

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f tmap-env.yml
```

## Input files

tMAP will work with any tRAX generated coverage file or the coverage file generated using this tool (using the --bedcoverage option)

### Metadata

You will need to provide the meta-data associated with the samples if you want to group samples based on specific experimental conditions. To do this you can provide a [tRAX samples file](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) whitespace/.csv/.tsv file (`--samples`) containing the sample names, sample groups, and path to fastq files. An example samples file is shown below:

```tsv
sample1 condition1 /path/to/sample1.fastq.gz
sample2 condition1 /path/to/sample2.fastq.gz
sample3 condition2 /path/to/sample3.fastq.gz
sample4 condition2 /path/to/sample4.fastq.gz
```

You will need to provide the meta-data associated with the comparisons you want to perform between experimental conditions. To do this you can provide a [tRAX pairs file](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) whitespace/.csv/.tsv file (`--exppairs`) containing the sample names, sample groups, and path to fastq files. An example samples file is shown below:

```tsv
condition1 condition2
condition1 condition3
condition2 condition3
```

## Usage

### Activating the environment

```bash
conda activate tmap
```

### Accessing the help information

```bash
python tMAP.py --h
```

### Determining potential sites of base modifications present in tRNAs

![alt text](https://github.com/Aimannin/TRMT1/tree/version-1.0/images/mismatch_coverage_example.jpg?raw=true)

The tMAP-getmismatches.py script can be used to identify potential sites of base modification present in tRNA sequencing data.

```bash
python tMAP-getmismatches.py \
--cov <tRAX_coverage_file> \
--o <output_file_prefix> \
[(optional) --org <organism> --alpha <mismatch_pseudocounts> --beta <reference_pseudocounts> --editfrac <minumum_edit_fraction> --minreads <minimum_read_coverage> --positions <list_of_sprinzl_positions> --multimapcutoff <percent_multimap_coverage> --exppairs <pairs.txt> --samples <samples.txt> --predict]
```
* `--cov` is the path to the tRAX coverage file you want to analyse modifications from; runname-coverage.txt.
* `--o`  is the prefix for the output files.
* (optional flags)
* `--positions` tRNA positions of interest. Can list multiple positions e.g. 9 26 32 34...
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--alpha` pseudocounts to add to mismatch counts for beta function; default=0.
* `--beta` pseudocounts to add to reference counts for beta function; default=0.
* `--editfrac` minimum fraction of mismatched reads for beta function; default=0.05.
* `--minreads` is the minumum number of read coverage of a base required to go into analysis; default=20.
* `--multimapcutoff` Percent of multimapping reads covering a site to mark as artefact; default=20.
* (optional flags; to determine differences in mismatch rates between all pairwise comparisons use these)
* `--exppairs` pairs.txt file from tRAX to perform pairwise comaprisons.
* `--samples` samples.txt file from tRAX, or custom made as previously described to group biological replicates.
* (optional flags (in development); to determine which base modification is present based on misincorporation signature)
* `--predict` (in development) predict which base modification is present based on misincorporation signature. Model built using tRNA sequencing (OTTR-seq) data from various eukaryotic model organisms.


### Getting isodecoder varation within an organism

The tMAP-align.py script can be used to generate an excel-friendly alignment of tRNA isodecoders, the base-level varation exising amongst isoacceptor groups (ex. all Arg-TCT sequence variation), and the pairwise sequence differences between tRNA isodecoders.

```bash
python tMAP-align.py \
--cov <tRAX_coverage_file> \
--alignments <dbname-trnaalign.stk> \
--o <output_file_prefix> \
[(optional) --minreads <minimum_read_coverage> --org <organism> --plot]
```

* `--cov` is the path to the tRAX coverage file you want to analyse modifications from.
* `--alignments` is the path to the stockholm file containing the mature tRNA alignments from tRAX makedb.py; dbname-trnaalign.stk.
* `--o`  is the prefix for the output files.
* (optional flags)
* `--minreads` is the minumum number of read coverage of a base required to go into analysis; default=20.
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
<!-- * `--plot` whether or not to generate a bar plot showing the Sprinzl positions containing instances of isodecoder variation. -->

### Getting base-resolution read coverage of user-defined features

The tMAP-bedcoverage.py script can be used to generate a tRAX-like coverage output file of user defined regions provided in bed format.

```bash
python tMAP-bedcoverage.py \
--fasta <genome.fa> \
--regions <regions_of_interest.bed> \
--alignments <sample1.bam> <sample2.bam> \
--o <output_file_prefix> \
[(optional) --minreads <minimum_read_coverage> --clipping <number_of_bases> --samples <samples.txt> --org <organism> --plot --threads <number_of_threads>]
```

* `--fasta` Genome fasta file; needs an index file (.fai), can be generated with samtools faidx fasta.fa.
* `--regions` bed file containing genomic region of interest.
* `--alignments` bam file(s) containing read alignments; needs index files (.fai), can be generated with samtools index alignments.bam. multiple bam files can be passed as follows; file1.bam file2.bam ...
* `--o`  is the prefix for the output files.
* (optional flags)
* `--minreads` is the minumum number of read coverage of a base required to go into analysis; default=20.
* `--clipping` umber of bases a read can extend beyond the feature start/end to be included in the analysis; default=0.
* `--samples` samples.txt file from tRAX, or custom made as previously described to group biological replicates.
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--plot` whether or not to generate coverage plots for each of the features.
* `--threads` number of threads to process the data on; Default=1.
