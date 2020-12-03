## Introduction
ERVsearch is a pipeline for identification of endogenous retrovirus like regions in a host genome, based on sequence similarity to known retroviruses.

ERVsearch screens for endogenous retrovirus (ERV) like regions in any FASTA file using the Exonerate algorithm (Slater and Birney, 2005, doi:[10.1186/1471-2105-6-31](https://doi.org/10.1186/1471-2105-6-31)). 

* In the **Screen** section, open reading frames (ORFs) resembling retroviral *gag*, *pol* and *env* genes are identified based on their level of similarity to a database of known complete or partial retroviral ORFs.
* In the **Classify** section, these ORFs are classified into groups based on a database of currently classified retroviruses and phylogenetic trees are built.
* In the **ERVRegions** section, regions with ORFs resembling more than one retroviral gene are identified.

This is a updated and expanded version of the pipeline used to identify ERVs in Brown and Tarlinton 2017 (doi: [10.1111/mam.12079](https://doi.org/10.1111/mam.12079)), Brown et al. 2014 (doi: [10.1128/JVI.00966-14](https://doi.org/10.1128/JVI.00966-14)), Brown et al. 2012 (doi: [j.virol.2012.07.010](https://doi.org/10.1016/j.virol.2012.07.010)) and Tarlinton et al. 2012 (doi: [10.1016/j.tvjl.2012.08.011](https://doi.org/10.1016/j.tvjl.2012.08.011)). The original version is available [here](https://github.com/ADAC-UoN/predict.genes.by.exonerate.pipeline) as a Perl pipeline and was written by Dr Richard Emes.


1. [Prerequisites](#prerequisites)<br>
2. [Installation](#installation)<br>
3. [Quick Start](#quick-start)<br>
4. [Pipeline Description](#pipeline-description)
5. [Input Files](#input-files)
6. [Usage](#usage)
7. [Parameters](parameters.html)
8. [Functions](functions.html)

## Prerequisites

The pipeline is currently available for Unix-based systems only.

The ERVsearch pipeline requires the following freely available software. All packages are available via pip and easy_install

Python 3.5+ with the following packages:
- [ruffus](https://pypi.python.org/pypi/ruffus)
- [numpy](https://pypi.python.org/pypi/numpy)
- [pandas](https://pypi.python.org/pypi/pandas)
- [ete3](https://pypi.python.org/pypi/ete3)
- [matplotlib](https://pypi.python.org/pypi/matplotlib)

The following commonly used software needs to be installed and in your $PATH
- [Samtools](https://sourceforge.net/projects/samtools)
- [Bedtools](https://github.com/arq5x/bedtools2)
- [Emboss](http://emboss.sourceforge.net)
- [Mafft](http://mafft.cbrc.jp/alignment/software)
- [FastTree](http://meta.microbesonline.org/fasttree)

The following software also needs to be installed
- [Exonerate](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
- [Usearch](http://www.drive5.com/usearch)

## Installation

The latest release of ERVsearch is available via pip3.

```
pip3 install ERVsearch
```

Alternatively, if you prefer to install directly,  the latest release can be downloaded from [Github](https://github.com/KatyBrown/ERVsearch/releases/latest).

The latest (beta) version can also be cloned from github

```
git clone https://github.com/KatyBrown/ERVsearch.git
```


No compliation is required, just add the ERVsearch directory to your path or use the full path to `ERVsearch/ERVsearch`. If installed using pip you should be able to call `ERVsearch` directly.


## Quick Start

After cloning the repository, the program can be used as is (with the above prerequisites installed).

1. Make a copy of the pipeline.ini file (ERVsearch/templates/pipeline.ini) in your working directory (the directory in which you would like to store the output).

2. Download a local copy of your genome (or other sequence) of interest as a single FASTA file.

3. Edit your copy of pipeline.ini to configure the pipeline for your computer:
* Add the path to the genome you want to screen (the fasta file in step 2) to the genome section
e.g. hg38.fa saved in /home/myname/genome/hg38.fa would require the following options:<br>

```
[genome]
file_=/home/myname/genome/hg38.fa
```

* Add the paths to ERVsearch, usearch and exonerate to the paths section
e.g.<br>
```
[paths]
path_to_ERVsearch=/home/myname/ERVsearch
path_to_usearch=/home/myname/usearch/usearch11.0.667_i86linux32
path_to_exonerate=/home/myname/exonerate/bin/exonerate
```

* Run the pipeline in your working directory as:
```
ERVsearch --target_tasks full -v5
```


## Pipeline Description
The pipeline is designed to identify regions resembling retroviral *gag*, *pol* and *env* genes in a genome (or other set of sequences) and to perform various analyses on these regions.

It is divided into three sections:

### Screen
[function documentation](functions.html#screen)<br>
* Screens the genome for ERV like regions by comparing the genome to a set of known retroviral ORFs using Exonerate.
* Confirms the Exonerate regions using UBLAST
* Finds and confirms ORFs within these regions
* Finds the most similar known retroviral ORF in the database to each of the newly identified ORFs<br>

### Classify
[function documentation](functions.html#classify)<br>
* Classifies the newly identified ORFs into groups based on the most similar known ORF
* Aligns the newly identified ORFs with reference sequences within these groups and builds a phylogenetic tree for each group.
* Finds clusters of newly identified ORFs within these trees
* Incorporates representative sequences from these clusters into a summary tree for each retroviral gene and genus  (based on classification into *gamma*, *beta*, *spuma*, *alpha*, *lenti*, *epsilon* and *delta* retroviruses as definied by the ICTV (https://talk.ictvonline.org/taxonomy).<br>

### ERVRegions
[function documentation](functions.html#ervregions)
* Identifies regions of the genome containing ORFs resembling more than one different retroviral gene within a certain distance<br>


All functions in all sections are described in detail in the [functions](functions.html) section.


## Input Files

### Required Input Files

1: *FASTA file to screen for ERVs*

The main input file is sequence file in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) containing DNA sequences from the genome which you wish to screen for ERV-like regions. This would usually be a reference or de novo assembled genome but can be any set of DNA sequences.

Reference genome sequences are available from [Ensembl](http://www.ensembl.org/info/about/species.html) and [UCSC](https://genome-euro.ucsc.edu/index.html) (amongst others). 

To be used as an input file, the reference genome needs to be contained in a single FASTA file.

For Ensembl genomes, this would usually be the `GENOMEID.dna.toplevel.fa.gz `file from the "download DNA sequence" page for the appropriate organism, substituting `GENOMEID` for the genome ID (e.g. GRCh38)

For UCSC genomes this would be `GENOMEID.fa.gz` from the bigZips directory for this organism on the FTP server, substituting `GENOMEID` for the genome ID (e.g. hg38)

It is possible to use a gzipped or zipped file, in which case the filename needs to end with .gz or .zip respectively.


2: `pipeline.ini` file

This file is a configuration file in [ini format](https://en.wikipedia.org/wiki/INI_file) containing the parmeters you wish to use.
This file needs to be in your working directory - the folder in which you wish to run ERVsearch.

A template `pipeline.ini` file should be used and edited - this file is available as `ERVsearch/templates/pipeline.ini`.

Options specified as `!?` are required, all others have a default value.

[Required parameters](parameters.html#required-parameters)

[Optional parameters](parameters.html#optional-parameters)


### Optional Input Files
#### Custom Databases
By default, ERVsearch will use the provided database of 774 ERV nucleotide sequences and corresponding amino acid sequences as a query against the provided genome. This database is designed to be representative of known retroviruses and to identify the majority of ERVs. However, a more specific custom database can also be provided and used for the initial screen.

To do this, the [database_use_custom_db](parameters.html#use-custom-db) parameter in the pipeline.ini can be set to `True`. The query sequences should be stored as FASTA files of amino acid sequences, with one file per retroviral gene. Only gag, pol and env genes are currently supported. Very short sequences (less than ~100 amino acids) should be avoided where possible.
The paths to these files are then specified in the database section of the pipeline.ini

e.g.
```
[database]
use_custom_db=True
gag=/home/katy/my_databases/gag_ervs.fasta
pol=/home/katy/my_databases/pol_ervs.fasta
env=/home/katy/my_databases/env_ervs.fasta
```

Currently, custom databases are only used for the initial Exonerate screen and UBLAST check, after this all classification based steps will use the default databases, as these sequences have been classified into subgroups for phylogenetic analysis.


#### Sequence List

`keep_chroms.txt`

A list of chromosome names to include. The names should the names in the fasta file cropped at the first space,  e.g. "NW_006711271.1 Panthera tigris altaica isolate TaeGuk unplaced genomic scaffold" should just be listed as NW_006711271.1. The names should be listed with one name per line, are case sensitive and need to be identical to those in the fasta file. This file needs to be named keep_chroms.txt and in the working directory.


## Usage
### Running the Pipeline
The pipeline is implemented using the pipeline development package [ruffus](http://www.ruffus.org.uk/), (Goodstadt 2010, doi:[10.1093/bioinformatics/btq524](https://doi.org/10.1093/bioinformatics/btq524)).

To run the full pipeline, the following command is used:
```
ERVsearch --target_tasks full
```

If the pipeline stops or fails at any point, it will restart from after the previous fully completed step (based on the presence of the output files from this step). If the output of an earlier step in the pipeline has a more recent timestamp than the output of a later step, the later step will be rerun.

Sections of the pipeline can be run as follows:

```
ERVsearch --target_tasks Screen
```
Screen with Exonerate and check the results with UBLAST, find ORFs and find the most similar known retroviral ORF.

```
ERVsearch --target_tasks Classify
```
Run the *Screen* steps and then sort the sequences into subgroups, build phylogenetic trees for these groups and a summary tree for each gene and genus.

```
ERVsearch --target_tasks ERVRegions
```
Run the *Screen* steps and then find regions with ORFs resembling more than one retroviral gene in close proximity.

```
ERVsearch --target_tasks full
```
Run all the above sections.

You can also run the pipeline up until any specific [function](functions.html) - any function name can be provided for the `target_tasks` parameter.

For example to run up until the end of the function `classifyWithExonerate`, use the following command.
```
ERVsearch --target_tasks classifyWithExonerate
```
All functions prior to this function will run if needed.


### Parallelisation
The pipeline is paralellised to run jobs simultaneously where possible. To do this, set the parameter`--jobs N` in the command line, where N is the number of CPUs available on your machine.
e.g.
```
ERVsearch --target_tasks full --jobs 8
```
This would run on 8 CPUs

If running on a high performance cluster, it is recommended to use a single node and set `--jobs` to the number of cores available on that node.

### Verbosity
Ruffus verbosity is set using the `-v` parameter from 1 to 10. The recommended setting for ERvsearch is -v 5, however this needs to be specified to override the ruffus default of 1.

e.g.
```
ERVsearch --target_tasks full -v 5
```
