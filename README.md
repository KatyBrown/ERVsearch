# ERVsearch

ERVsearch is a pipeline for identification of endogenous retrovirus like regions in a host genome, based on sequence similarity to known retroviruses.


## Introduction

ERVsearch screens for endogenous retrovirus (ERV) like regions in any FASTA file using the Exonerate algorithm (Slater and Birney, 2005, doi 10.1186/1471-2105-6-31). 

* In the **Screen** section, open reading frames (ORFs) resembling retroviral *gag*, *pol* and *env* genes are identified based on their level of similarity to a database of known complete or partial retroviral ORFs.
* In the **Classify** section, these ORFs are classified into groups based on a database of currently classified retroviruses.
* In the **ERVRegions** section, regions with ORFs resembling more than one retroviral gene are identified.
* In the *Summarise** section various summary plots are produced.

## Prerequisites

The pipeline is currently available for Unix-based systems only.

The ERVsearch pipeline requires the following freely available software. All packages are available via pip and easy_install

Python 3.5+ with the following packages:
- ruffus - https://pypi.python.org/pypi/ruffus
- numpy - https://pypi.python.org/pypi/numpy
- pandas - https://pypi.python.org/pypi/pandas
- ete3 - https://pypi.python.org/pypi/ete3
- matplotlib - https://pypi.python.org/pypi/matplotlib

The following commonly used software needs to be installed and in your $PATH
- Samtools: https://sourceforge.net/projects/samtools/files/
- Bedtools: https://github.com/arq5x/bedtools2
- Emboss: http://emboss.sourceforge.net/download/#Stable/
- Mafft: http://mafft.cbrc.jp/alignment/software/linux.html
- FastTree: http://meta.microbesonline.org/fasttree/#Install

The following software also needs to be installed
- Exonerate: http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
- Usearch: http://www.drive5.com/usearch/download.html

Installation
--------------
The latest release of ERVsearch is available via pip3.
`pip3 install ERVsearch`

Alternatively, if you prefer to install directly,  the latest release can be downloaded from Github.
https://github.com/KatyBrown/ERVsearch/releases/latest

The latest (beta) version can also be cloned from github
`git clone https://github.com/KatyBrown/ERVsearch.git`

No compliation is required, just add the ERVsearch directory to your path or use the full path to ERVsearch/ERVsearch. If installed using pip you should be able to call `ERVsearch` directly.


## Quick Start

After cloning the repository, the program can be used as is (with the above prerequisites installed).

1. Make a copy of the pipeline.ini file (templates/pipeline.ini) in your working directory (the directory in which would would like to store the output).

2. Download a local copy of your genome (or other sequence) of interest as a single FASTA file.

3. Edit your copy of pipeline.ini to configure the pipeline for your computer:
	* Add the path to the genome you want to screen (the fasta file in step 2) to the genome section
	 	 e.g. hg38.fa saved in /home/myname/genome/hg38.fa would require the following options:
	 	`[genome]`
	 	`file_=/home/myname/genome/hg38.fa`

	* Add the paths to ERVsearch, usearch and exonerate to the paths section
		e.g.
		`[paths]`
		`path_to_ERVsearch=/home/myname/ERVsearch`
		`path_to_usearch=/home/myname/usearch/usearch11.0.667_i86linux32`
		`path_to_exonerate=/home/myname/exonerate/bin/exonerate`
		
	* Run the pipeline in your working directory as:
		`ERVsearch --target_tasks full -v5` 


##Input Files

###Required Input Files

1: *FASTA file to screen for ERVs*

The main input file is sequence file in FASTA format (https://en.wikipedia.org/wiki/FASTA_format) containing DNA sequences from the genome of interest which you wish to screen for ERV-like regions. This would usually be a reference or de novo assembled genome but can be any set of DNA sequences.

Reference genome sequences are available from http://www.ensembl.org/info/about/species.html and https://genome-euro.ucsc.edu/index.html (amongst others). 

To be used as an input file, the reference genome needs to be contained in a single FASTA file.

For Ensembl genomes, this would usually be the` GENOMEID.dna.toplevel.fa.gz `file from the "download DNA sequence" page for the appropriate organism, substituting `GENOMEID` for the genome ID (e.g. GRCh38)

For UCSC genomes this would be `GENOMEID.fa.gz` from the bigZips directory for this organism on the FTP server, substituting `GENOMEID` for the genome ID (e.g. hg38)

It is possible to use a gzipped or zipped file, in which case the filename needs to end with .gz or .zip respectively.


2: `pipeline.ini` file

This file is a configuration file in `ini` format (https://en.wikipedia.org/wiki/INI_file) containing the parmeters you wish to use.
This file needs to be in your working directory - the folder in which you wish to run ERVsearch.

A template `pipeline.ini` file should be used and edited - this file is available as templates/pipeline.ini.

Options specified as `!?` are required, all others have a default value.

Required parameters are as follows:

|Section|Parameter|Description|Example|
|-|-|-|-|
|`genome`|`file`|path to genome to screen|`/home/katy/ERVsearch_screen/mygenome.fasta`|
|`paths`|`path_to_ERVsearch`| path to the main ERVsearch folder on your computer|`/home/katy/ERVsearch`|
|`paths`|`path_to_exonerate`|path to the exonerate executable on your computer|`/home/katy/Exonerate2.4.0/bin/exonerate
|`paths`|`path_to_usearch`|path to the USEARCH executable on your computer| `/home/katy/USEARCH/usearch`|

Optional parameters are listed in the **Parameters** section below.


###Optional Input Files

`keep_chroms.txt`

A list of chromosome names to include. The names should the names in the fasta file cropped at the first space,  e.g. "NW_006711271.1 Panthera tigris altaica isolate TaeGuk unplaced genomic scaffold" should just be listed as NW_006711271.1. The names should be listed with one name per line, are case sensitive and need to be identical to those in the fasta file. This file needs to be named keep_chroms.txt and in the working directory.


## Usage

### Screen
#### Parameters
#### Outputs
### Classify
#### Parameters
#### Outputs
### ERVRegions
#### Parameters
#### Outputs
### Summarise
#### Parameters
#### Outputs




## Pipeline Functions
### Screen
### Classify
### ERVRegions
### Summarise

## Parameters
