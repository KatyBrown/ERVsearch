#!/usr/bin/env python3
'''
Functions to identify regions with multiple ERV-like ORFs
'''
import pandas as pd
import HelperFunctions
import ORFs


def getDicts(genes, mergecols):
    '''
    mergecols is the information from the grouped.dir/GENE_groups.tsv which
    is relevant in the ERV_regions output table - the most similar reference
    sequence ("match" column), group and genus.

    Returns a dictionary of dictionaries with a key for each gene, for each
    column in mergecols geneD[gene][column] is a dictionary where
    keys are ERV IDs and values are the value of the column in the
    grouped.dir/GENE_groups.tsv table.
    '''
    geneD = dict()
    for gene in genes:
        genetab = pd.read_csv("grouped.dir/%s_groups.tsv" % gene, sep="\t")
        geneD[gene] = dict()
        for col in mergecols:
            geneD[gene][col] = dict(zip(genetab['ID'], genetab[col]))
    return (geneD)


def makeColumns(genes):
    '''
    Get the column names for the various tables
    '''
    # standard bed file columns
    bed_cols = HelperFunctions.getBedColumns()

    # rename the "name" column as "orig_name" because it will be merged
    # with a table where "name" exists but is different
    bed_cols[3] = 'orig_name'

    # columns from the grouped.dir/GENE_groups.tsv table which have relevant
    # information
    # match is the most similar known retrovirus, group is the assigned
    # retroviral group, genus is the retroviral genus
    merge_cols = ['match', 'group', 'genus']

    # temporary for storing columns to add to keep_cols, final_cols and
    # all_cols
    new_cols = ['name', 'ID', 'start', 'end', 'strand'] + merge_cols

    # all columns in the final table
    all_cols = ['name']

    # columns to keep in the output table

    keep_cols = ['name', 'chrom', 'start', 'end', 'strand', 'genus']
    for gene in genes:
        for col in new_cols:
            all_cols.append("%s_%s" % (gene, col))
            keep_cols.append("%s_%s" % (gene, col))
            if col == 'start' or col == "end":
                keep_cols.append("%s_relative_%s" % (gene, col))
    all_cols.append('orig_name')
    keep_cols.append('orig_name')
    return (bed_cols, merge_cols, all_cols,  keep_cols)


def checkRegion(genes, region):
    '''
    Check which genes are present in a region based on a string containing
    the concatenated "name" column from the bedtools merged output - the
    names of the merged regions delimited by ",".
    e.g. pol247_NW_006711412.1-501724-502357(-),gag154_NW_006711412.1-502959-503325(-)
    contains pol and gag.
    '''
    region_genes = []
    for gene in genes:
        if gene in region:
            region_genes.append(gene)
    return (region_genes)


def processRegion(region, gene, merge_cols, geneD):
    '''
    Take all the instances of a particular gene found in a region and process
    them into columns for the overall regions table.
    Often there will be more than one ORF associated with a gene, especially
    the longer genes (usually pol) - this function merges these into a
    single output.
    Generates a list of 9 results:
        * new_name: the name (ID_chrom_start-end(strand)) of the gene region
        defined here
        * new_ID: the IDs of the instances of this gene in this region
        * new_start: the minimum start position of all instances of this
        gene in this region
        * new_end: the maximum end position of all instances of this gene in
        this region
        * new_strand: the strand or strands of this gene in this region
        * new_match: the reference sequences matched to these instances of
        this gene using Exonerate
        * new_group: the group these instances of this gene were assigned to
        * new_genus: the genera these instances of this gene were assigned to

    This is carried out using the the concatenated "name" column from the
    bedtools merged output - the names of the merged regions delimited by ",".
    '''
    # split on the delimiter
    R = region.split(",")
    # is this gene there at all?
    this_gene = []
    for seg in R:
        if gene in seg:
            this_gene.append(seg)
    # if this gene is present in this region
    if len(this_gene) != 0:
        # store the IDs, starts, ends, strands
        IDs, starts, ends, strands = [], [], [], []
        # dictionary to store the group information for this gene
        D = {c: [] for c in merge_cols}
        for ID in this_gene:
            # get the ORF ID, chromosome, start, end and strand
            namD = ORFs.splitNam(ID)
            starts.append(namD['start'])
            ends.append(namD['end'])
            IDs.append(namD['ID'])
            strands.append(namD['strand'])
            for col in merge_cols:
                D[col].append(geneD[gene][col][namD['ID']])
        # find the minimum start positoin
        new_start = min(starts)
        # find the maximum end position
        new_end = max(ends)
        # concatenate the IDs with | as a delimiter
        new_ID = "|".join(IDs)
        # concatenate the strands (unique only)
        new_strand = "".join(uniqueList(strands))

        # Make a new ID based on this info
        new_name = "%s_%s-%s-%s(%s)" % (new_ID, namD['chrom'],
                                        new_start, new_end,
                                        new_strand)

        # combine this to make part of the table row
        new_row = [new_name, new_ID, new_start, new_end, new_strand]

        # add the group information to the row
        for col in merge_cols:
            new_row.append("|".join(uniqueList(D[col])))
    else:
        # if this gene isn't found in this region, use NA as a placeholder
        # so the final output table has the same columns in every row.
        # float('nan') is used for numeric columns
        new_row = ["NA", "NA", float('nan'), float('nan'),
                   "NA", "NA", "NA", "NA"]
    return (new_row)


def uniqueList(x):
    '''
    Take a list and return the unique values (as a list)
    '''
    return (list(set(x)))


def fixCol(genes, results, col):
    '''
    Combine the columns labeled gene_col for each gene without duplicates or
    NAs
    '''
    G = results['%s_%s' % (genes[0], col)] + "|"
    for gene in genes[1:]:
        G += results['%s_%s' % (gene, col)] + "|"
    new_G = []
    for g in G:
        g_split = uniqueList(g.split("|"))
        g_merged = "|".join(g_split)
        g_final = g_merged.replace("NA", "").replace("||", "").strip("|")
        new_G.append(g_final)
    return (new_G)


def cleanTable(results, bed, keep_cols, genes):
    '''
    Tidy up the final output table:
        * Merge with the original bed file to get the start and end
          co-ordinates of the region
        * Sort by chromosome then start position
        * Merge strand and genus columns to make a combined strand and genus
          for all genes
        * Find the start and end positions of the genes relative to the
          start of the region (rather than the chromosome positions)
    for all genes, remove ununsed columns.

    '''
    results = bed.merge(results)
    results = results.sort_values(['name', 'chrom', 'start'])
    results['strand'] = fixCol(genes, results, 'strand')
    results['genus'] = fixCol(genes, results, 'genus')
    for gene in genes:
        results['%s_relative_start' % gene] = results[
            '%s_start' % gene] - results['start']
        results['%s_relative_end' % gene] = results[
            "%s_end" % gene] - results['start']
    results = results[keep_cols]
    results = results.fillna("NA")
    return (results)


def getRegions(infile, genes, log):
    '''
    Takes a merged bed file consisting of regions of the genome identified
    as having more than one ERV-like ORF, finds the regions within this file
    which contain more than one different gene (e.g. gag and pol instead of
    two gag ORFs) and outputs a formatted table of information about these
    regions.

    The output table will usually have 37 columns:
        * name - the final ID of the ERV region - the genes found plus an
          integer
          e.g. gag_pol_12
        * chrom - chromosome
        * start - start position of the ERV region
        * end - end position of the ERV region
        * strand - strand of the ERv region
        * genus - genus of the ERV region, can be multiple genera delimted by
          | if different genes had different genera
        * for each gene screened for (usually gag, pol and env)
              * GENE_name - the names of the ORFs for this gene in this region
              * GENE_ID - the original IDs of the ORFs for this gene in this
                region
              * GENE_start - the start position of this gene in this region
                (genome co-ordinates)
              * GENE_relative_start - the start position of this gene in this
                region (relative to the start of the region)
              * GENE_end - the end position of this gene in this region
                (genome co-ordinates)
              * GENE_relative_end - the end position of this gene in this
                region (relative to the end of the region)
              * GENE_strand - the strand for this gene in this region
              * GENE_match - the closest reference retrovirus to this gene
                in this region
              * GENE_group - the group of the closest reference retrovirus to
                this gene in this region
              * GENE_genus - the genus of the closest reference retrovirus to
                this gene in this region
        * orig_name - the name of the region in the input table
    '''
    log.info("Retrieving column names for %s" % infile)
    # get the column names for all the tables
    bed_cols, merge_cols, all_cols, keep_cols = makeColumns(genes)

    # dictionaries containing the most similar known retrovirus, group and
    # genus of the ERV regions in the input table
    geneD = getDicts(genes, merge_cols)

    # read the bed file into a table
    bed = pd.read_csv(infile, sep="\t", header=None)
    # add column names
    bed.columns = bed_cols

    # empty variable to store the output table
    rows = []

    # this is used to keep count of how many regions of each type have been
    # found so they can be named
    regionD = dict()
    log.info("Annotating multi-gene regions in %s" % infile)
    for region in bed['orig_name']:
        # find out which genes are in the region based on the names of the
        # ORFs in the region
        region_genes = checkRegion(genes, region)
        # if there is more than one different gene in the region
        if len(region_genes) > 1:
            # the regions are named with the concatenated names of the
            # genes which are present as a prefix e.g. gag_pol
            region_gene_ID = "_".join(region_genes)
            regionD.setdefault(region_gene_ID, 0)
            # keep count of the number of regions with the same prefix
            regionD[region_gene_ID] += 1
            # assign an ID - the prefix plus the next number
            region_ID = "%s_%i" % (region_gene_ID, regionD[region_gene_ID])
            # the first column in the output should be the region ID
            row = [region_ID]
            for gene in genes:
                # Get the details of this gene in this region to add to the
                # row about this region
                gene_results = processRegion(region, gene, merge_cols, geneD)
                row += gene_results
            row.append(region)
            rows.append(row)
    # make a data frame of all the rows
    results = pd.DataFrame(rows, columns=all_cols)
    log.info("Identified %i multi-gene regions in %s" % (len(results), infile))
    # tidy up the data frame
    clean_results = cleanTable(results, bed, keep_cols, genes)
    return (clean_results)
