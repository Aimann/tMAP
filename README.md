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

### getmismatch -- Determining potential sites of base modifications present in tRNAs

The getmismatch module can be used to identify potential sites of base modification present in tRNA sequencing data. Identification of potentially modified bases was performed by building on the Bayesian method for genomic variant calling ([Li et al., 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2577856/); [Washburn et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959997/); [Deffit et al., 2017](https://elifesciences.org/articles/28625)) with a custom-designed prior on the % mismatches to account for sequencing errors. Pairwise comparisons of mismatch rates can be determined by providing a [tRAX-like](http://trna.ucsc.edu/tRAX/#step-3-analyze-sequencing-data-for-gene-expression) samples file as well as a pairs file. P-values are computed for each site using the T-test, and adjusted p-values were determined using Benjamini/Hochberg multiple-test correction.

```bash
python tMAP.py getmismatch\
--cov <tRAX_coverage_file> \
--o <output_file_prefix> \
[(optional) --positions <list_of_sprinzl_positions> --exppairs <pairs.txt> --samples <samples.txt> --org <organism> --alpha <mismatch_pseudocounts> --beta <reference_pseudocounts> --editfrac <minumum_edit_fraction> --minreads <minimum_read_coverage> --multimapcutoff <percent_multimap_coverage> --predict]
```
* `--cov` is the path to the tRAX coverage file you want to analyze modifications from; runname-coverage.txt.
* `--o`  is the prefix for the output files.

(optional flags)

* `--positions` tRNA positions of interest. Can list multiple positions e.g. 9 26 32 34...
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--alpha` pseudocounts to add to mismatch counts for beta function; default=0.
* `--beta` pseudocounts to add to reference counts for beta function; default=0.
* `--editfrac` minimum fraction of mismatched reads for beta function; default=0.05.
* `--minreads` is the minimum number of read coverage of a base required to go into the analysis; default=5.
* `--multimapcutoff` Percent of multimapping reads covering a site to mark as artifact; default=20.
  
(optional flags; to determine differences in mismatch rates between all pairwise comparisons, use these)
  
* `--exppairs` pairs.txt file from tRAX to perform pairwise comaprisons.
* `--samples` samples.txt file from tRAX, or custom made as previously described to group biological replicates.
* (optional flags (in development); to determine which base modification is present based on misincorporation signature)
* `--predict` (in development) predict which base modification is present based on misincorporation signature. Model built using tRNA sequencing (OTTR-seq) data from various eukaryotic model organisms.


### align -- Getting isodecoder variation within an organism

The align module can be used to generate an excel-friendly alignment of tRNA isodecoders, the base-level variation existing amongst isoacceptor groups (e.g., all Arg-TCT sequence variation), and the pairwise sequence differences between tRNA isodecoders (e.g. Arg-TCT-1 vs Arg-TCT-2)

```bash
python tMAP.py align \
--stk <dbname-trnaalign.stk> \
--o <output_file_prefix> \
[(optional) --org <organism>]
```

* `--alignments` is the path to the Stockholm file containing the mature tRNA alignments from tRAX makedb.py; dbname-trnaalign.stk.
* `--o`  is the prefix for the output files.
  
(optional flags)

* `--org` organism of interest (euk, bact, arch, mito); default='euk'.

### bedregion -- Getting base-resolution read coverage of user-defined features

The bedregion module can be used to generate a tRAX-like coverage output file of user-defined regions provided in bed format. Useful for looking at misincorporation data in any non-tRNA features like rRNA, miRNAs, snoRNAs, etc.

```bash
python tMAP.py bedregion \
--fasta <genome.fa> \
--regions <regions_of_interest.bed> \
--alignments <sample1.bam> <sample2.bam> \
--o <output_file_prefix> \
[(optional) --minreads <minimum_read_coverage> --clipping <number_of_bases> --samples <samples.txt> --org <organism> --plot --threads <number_of_threads>]
```

* `--fasta` Genome fasta file; needs an index file (.fai), can be generated with samtools faidx fasta.fa.
* `--regions` bed file containing genomic regions of interest. If using BED12 format with multiple exons also use the `--split` flag
* `--alignments` bam file(s) containing read alignments; needs index files (.fai), can be generated with samtools index alignments.bam. Multiple bam files can be passed one after the other as follows; file1.bam file2.bam ...
* `--o`  is the prefix for the output files.
  
(optional flags)

* * `--threads` the number of threads to process the data on; Default=1.
* `--clipping` the number of bases a read can extend beyond the feature start/end to be included in the analysis; default=2.
* `--samples` samples.txt file from tRAX, or custom made as previously described to group biological replicates.
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--plot` whether or not to generate coverage plots for each of the features.
* `--split` whether or not using BED12 format; for transcripts with multiple exons will concatenate them together.
