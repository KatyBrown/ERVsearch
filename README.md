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
After cloning the repository, the program can be used as is (with the above prerequisites installed).

1. Make a copy of the pipeline.ini file (src/pipeline_ERVs/pipeline.ini) in your working directory (the directory in which would would like to store the output).

2. Download a local copy of your genome of interest as a single fasta file.  Check if the genome is assembled into chromosomes, scaffolds or contigs.

3. Edit this copy to configure the pipeline for your computer:

    Change the working directory to the absolute path to the working directory.

    Add the name of your genome and the absolute path to the directory where you have saved it

      e.g. hg19.fa saved in /home/myname/genome/hg19.fa would require the following options:
      
        - genome=hg19
        - genome_directory=/home/myname/genome
    
    If your genome is assembled into chromosomes, has_chroms should be 1, otherwise it should be 0.

    Change the paths to the ERV input files to the directories in which you have cloned the repository.

      e.g if you cloned to /home/ERVsearch then paths would be
      
        - path_to_refs=/home/ERVsearch/ERV_db/all_ERVs.fasta
        - path_to_phyloseqs=/home/ERVsearch/phylogenies
        - sequencedir=/home/ERVsearch/ERV_db

    Change the paths to the required software.

4. Run the pipeline in the working directory as:

    -python /path_to_ERVsearch/ERVsearch/src/pipeline_ERVs.py -v


Details
-------
Input
-----
Output
------
