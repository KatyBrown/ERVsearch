# ERVsearch

Full documentation is available via [ReadTheDocs](http://ervsearch.readthedocs.io).

ERVsearch is a pipeline for identification of endogenous retrovirus like regions in a host genome, based on sequence similarity to known retroviruses.

ERVsearch screens for endogenous retrovirus (ERV) like regions in any FASTA file using the Exonerate algorithm (Slater and Birney, 2005, doi:[10.1186/1471-2105-6-31](https://doi.org/10.1186/1471-2105-6-31)). 

* In the **Screen** section, open reading frames (ORFs) resembling retroviral *gag*, *pol* and *env* genes are identified based on their level of similarity to a database of known complete or partial retroviral ORFs.
* In the **Classify** section, these ORFs are classified into groups based on a database of currently classified retroviruses and phylogenetic trees are built.
* In the **ERVRegions** section, regions with ORFs resembling more than one retroviral gene are identified.

This is a updated version of the pipeline used to identify ERVs in Brown et al. 2014 (doi: [10.1128/JVI.00966-14](https://doi.org/10.1128/JVI.00966-14)) and Brown and Tarlinton 2016 (doi: [10.1111/mam.12079](https://doi.org/10.1111/mam.12079)). The original version is available [here](https://github.com/ADAC-UoN/predict.genes.by.exonerate.pipeline) as a Perl pipeline and was written by Dr Richard Emes.