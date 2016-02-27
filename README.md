# ERVsearch

ERVsearch is a ruffus pipeline for identification of endogenous retrovirus like regions in a host genome.

Introduction
------------
ERVsearch provides a comprehensive screen for endogenous retrovirus (ERV) like regions in any genome.  The genome needs to be available as a single FASTA file.
Reference genome sequences are available from http://www.ensembl.org/info/about/species.html and https://genome-euro.ucsc.edu/index.html. De novo assemblies can also be used.

ERV search will take this genome and output:
A table showing regions of the genome related to each of the major retroviral genes (gag, pol and env).  This table shows:
- ID, chromosome, start position, end position and length of the ERV region
- The name and source of the most similar previously known retrovirus sequence and its similarity to the region
- The start position, end position and length of the longest open reading frame (ORF) in the region

Summary plots for each gene, consisting of:
- A bar plot showing the number of genes per chromosome
- A histogram of the ORF length
- A pie chart showing distribution between retroviral genera
- A bar plot showing the frequencies of regions related to different known retroviruses

Prerequisites
-------------
The pipeline is currently available for Unix systems only.

The ERVsearch pipeline requires the following freely available software.

Python 2 with the following packages:
- ruffus - https://pypi.python.org/pypi/ruffus
- numpy - https://pypi.python.org/pypi/numpy
- pandas - https://pypi.python.org/pypi/pandas
- ete2 - https://pypi.python.org/pypi/ete2
- matplotlib - https://pypi.python.org/pypi/matplotlib

All packages are available via pip and easy_install


- Exonerate: http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
- Bedtools: https://github.com/arq5x/bedtools2
- Samtools: https://sourceforge.net/projects/samtools/files/
- Usearch: http://www.drive5.com/usearch/download.html
- Emboss: http://emboss.sourceforge.net/download/#Stable/
- Mafft: http://mafft.cbrc.jp/alignment/software/linux.html
- FastTree: http://meta.microbesonline.org/fasttree/#Install


Quick Start
-----------
Details
-------
Input
-----
Output
------
