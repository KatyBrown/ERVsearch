#!/usr/bin/env python3
'''
Functions to identify regions with multiple ERV-like ORFs
'''
import pandas as pd
import HelperFunctions
import ORFs
import matplotlib
import matplotlib.pyplot as plt
import itertools
import numpy as np


def plotERVRegions(table, genes, plotparams, log):
    dpi = int(plotparams['dpi'])
    for ind in table.index.values:
        f = plt.figure(figsize=(10, 3), dpi=dpi)
        a = f.add_subplot(111)
        row = table.loc[ind]
        a.hlines(1.5, row['start'], row['end'], zorder=0)
        length = row['end'] - row['start']
        a.set_xlim((row['start'] - (length * 0.01)),
                   (row['end'] + (length * 0.01)))
        a.set_ylim(0, len(genes) + 1)
        a.ticklabel_format(useOffset=False, style='plain')
        ticks = []
        for j, gene in enumerate(genes):
            if isinstance(row['%s_ID' % gene], str):
                gene_start = row['%s_start' % gene]
                gene_end = row['%s_end' % gene]
                gene_length = gene_end - gene_start
                gene_colour = plotparams['%s_colour' % gene]
                gene_name = row['%s_ID' % gene]
                pos = gene_end - (gene_length / 2)
                a.add_patch(matplotlib.patches.Rectangle((gene_start, 1+j),
                                                         gene_length,
                                                         1, color=gene_colour))
                a.text(pos, 1.5+j, gene_name, ha='center', va='center')
                ticks.append(int(gene_start))
                ticks.append(int(gene_end))

        a.set_xticks(ticks)
        a.set_xticklabels(ticks, rotation='vertical')
        a.set_title("%s, %s nts (%s)" % (row['name'], length, row['strand']))
        a.set_yticks([])
        f.savefig("ERV_region_plots.dir/%s.%s" % (row['name'],
                                                  plotparams['format']),
                  dpi=dpi, bbox_inches='tight')
        plt.close()


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


def overlap(start1, end1, start2, end2):
    """
    Check if the two ranges overlap.
    """
    return end1 >= start2 and end2 >= start1


def checkOverlaps(row, genes, max_overlap):
    '''
    Excludes genes where the overlapping region is a proportion > max_overlap
    of either gene - for some of the reference sequences there's not a clear
    delineation between especially between gag and pol.
    '''
    max_overlap = float(max_overlap)
    keep_genes = []
    # relative positions of start and end in the columns for each gene
    start_pos = 2
    end_pos = 3
    for i, j in itertools.combinations(np.arange(len(genes)), 2):
        section1 = [((i * 8) + 1), ((i+1) * 8) + 1]
        section2 = [((j * 8) + 1), ((j+1) * 8) + 1]

        gene1_start = row[section1[0] + start_pos]
        gene1_end = row[section1[0] + end_pos]
        gene2_start = row[section2[0] + start_pos]
        gene2_end = row[section2[0] + end_pos]

        if not np.isnan(gene1_start) and not np.isnan(gene2_start):
            if gene1_start < gene2_start:
                overlap = gene2_start - gene1_end
            else:
                overlap = gene1_start - gene2_end
            if overlap <= 0:
                p1 = np.abs(overlap) / (gene1_end - gene1_start)
                p2 = np.abs(overlap) / (gene2_end - gene2_start)
                if p1 > p2 and p1 > max_overlap:
                    row_left = row[:section1[0]]
                    row_right = row[section1[1]:]
                    row = row_left + ['NA', 'NA', float('nan'), float('nan'),
                                      "NA", "NA", "NA", "NA"] + row_right
                    keep_genes.append(genes[j])
                elif p2 >= p1 and p2 > max_overlap:
                    row_left = row[:section2[0]]
                    row_right = row[section2[1]:]
                    row = row_left + ['NA', 'NA', float('nan'), float('nan'),
                                      "NA", "NA", "NA", "NA"] + row_right
                    keep_genes.append(genes[i])
                else:
                    keep_genes.append(genes[i])
                    keep_genes.append(genes[j])
            else:
                keep_genes.append(genes[i])
                keep_genes.append(genes[j])
    keep_genes = set(keep_genes)
    return (row, keep_genes)


def getRegions(infile, genes, max_overlap, log):
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
            new_row, new_genes = checkOverlaps(row, genes, max_overlap)
            if len(new_genes) > 1:
                if len(new_genes) == len(region_genes):
                    row.append(region)
                    rows.append(new_row)
                else:
                    region_genes2 = []
                    for g in region_genes:
                        if g in new_genes:
                            region_genes2.append(g)
                    new_region_gene_ID = "_".join(region_genes2)
                    regionD[region_gene_ID] -= 1
                    regionD.setdefault(new_region_gene_ID, 0)
                    regionD[new_region_gene_ID] += 1
                    new_region_ID = "%s_%i" % (new_region_gene_ID,
                                               regionD[new_region_gene_ID])
                    row[0] = new_region_ID
                    print ("help")
            else:
                regionD[region_gene_ID] -= 1
    # make a data frame of all the rows
    results = pd.DataFrame(rows, columns=all_cols)
    log.info("Identified %i multi-gene regions in %s" % (len(results), infile))
    # tidy up the data frame
    clean_results = cleanTable(results, bed, keep_cols, genes)
    return (clean_results)
