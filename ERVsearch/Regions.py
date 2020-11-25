#!/usr/bin/env python3
'''
Functions to identify regions with multiple ERV-like ORFs
'''
import pandas as pd
import HelperFunctions
import ORFs


def getDicts(genes, mergecols):
    '''
    Make a dictionary of dictionaries with a key for each gene.
    For each column in mergecols geneD[gene][column] is a dictionary where
    keys are ERV IDs and values are the value of the column in the
    grouped.dir/gene_groups.tsv table.
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
    bed_cols = HelperFunctions.getBedColumns()
    bed_cols[3] = 'orig_name'
    genes = ['gag', 'pol', 'env']
    merge_cols = ['match', 'group', 'genus']
    new_cols = ['name', 'ID', 'start', 'end', 'strand'] + merge_cols
    all_cols = ['name']
    final_cols = []
    keep_cols = ['name', 'chrom', 'start', 'end', 'strand', 'genus']
    for gene in genes:
        for col in new_cols:
            all_cols.append("%s_%s" % (gene, col))
            final_cols.append("%s_%s" % (gene, col))
            keep_cols.append("%s_%s" % (gene, col))
            if col == 'start' or col == "end":
                keep_cols.append("%s_relative_%s" % (gene, col))
    all_cols.append('orig_name')
    keep_cols.append('orig_name')
    return (bed_cols, merge_cols, all_cols, final_cols, keep_cols)


def checkRegion(genes, region):
    '''
    Check which genes are present in a region.
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
    Generates a list of 9 results:
        new_name: the name (ID_chrom_start-end(strand)) of the gene region
        defined here
        new_ID: the IDs of the instances of this gene in this region
        new_start: the minimum start position of all instances of this
        gene in this region
        new_end: the maximum end position of all instances of this gene in
        this region
        new_strand: the strand or strands of this gene in this region
        new_match: the reference sequences matched to these instances of
        this gene using Exonerate
        new_group: the group these instances of this gene were assigned to
        new_genus: the genera these instance of this gene were assigned to
    '''
    R = region.split(",")
    this_gene = []
    for seg in R:
        if gene in seg:
            this_gene.append(seg)
    if len(this_gene) != 0:
        IDs, starts, ends, strands = [], [], [], []
        D = {c: [] for c in merge_cols}
        for ID in this_gene:
            namD = ORFs.splitNam(ID)
            starts.append(namD['start'])
            ends.append(namD['end'])
            IDs.append(namD['ID'])
            strands.append(namD['strand'])
            for col in merge_cols:
                D[col].append(geneD[gene][col][namD['ID']])
        new_start = min(starts)
        new_end = max(ends)
        new_ID = "|".join(IDs)
        new_strand = "".join(uniqueList(strands))

        new_name = "%s_%s-%s-%s(%s)" % (new_ID, namD['chrom'],
                                        new_start, new_end,
                                        new_strand)
        new_row = [new_name, new_ID, new_start, new_end, new_strand]

        for col in merge_cols:
            new_row.append("|".join(uniqueList(D[col])))
    else:
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


def cleanTable(results, bed, final_cols, keep_cols, genes):
    '''
    Tidy up the final output table - sort, merge strand and genus columns
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


def getRegions(infile, genes):
    bed_cols, merge_cols, all_cols, final_cols, keep_cols = makeColumns(genes)
    geneD = getDicts(genes, merge_cols)
    bed = pd.read_csv(infile, sep="\t", header=None)
    bed.columns = bed_cols
    rows = []
    regionD = dict()

    for region in bed['orig_name']:
        region_genes = checkRegion(genes, region)
        # if there is more than one different gene in the region
        if len(region_genes) > 1:
            region_gene_ID = "_".join(region_genes)
            regionD.setdefault(region_gene_ID, 0)
            regionD[region_gene_ID] += 1
            region_ID = "%s_%i" % (region_gene_ID, regionD[region_gene_ID])
            row = [region_ID]
            for gene in genes:
                gene_results = processRegion(region, gene, merge_cols, geneD)
                row += gene_results
            row.append(region)
            rows.append(row)
    results = pd.DataFrame(rows, columns=all_cols)
    clean_results = cleanTable(results, bed, final_cols, keep_cols, genes)
    return (clean_results)
