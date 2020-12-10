#!/usr/bin/env python3
"""
Functions for generating summary tables and plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import HelperFunctions
import Fasta
import ORFs


def getInfiles(allinfiles, gene, ind):
    '''
    Get the input files corresponding to a particular index and gene
    '''
    if isinstance(allinfiles[0], list):
        inf = np.array([i[ind] for i in allinfiles])
    else:
        inf = np.array(allinfiles)
    inf_stems = np.array([os.path.basename(x) for x in inf])
    infile = inf[np.char.startswith(inf_stems, gene)][0]
    return (infile)


def getOutfiles(D, plot_format):
    outfiles = []
    for section in D:
        o = []
        for suffix in D[section]:
            if suffix.endswith(".txt"):
                o.append("summary_tables.dir/%s_%s" % (section,
                                                       suffix)),
            else:
                o.append("summary_plots.dir/%s_%s.%s" % (section,
                                                         suffix,
                                                         plot_format))
        outfiles.append(o)
    return (outfiles)


def allScreenOutfiles(plot_format):
    D = {'exonerate_initial': ['summary.txt',
                               'lengths', 'scores', 'strands',
                               'by_sequence', 'counts_per_gene'],
         'ublast_hits': ['initial_summary.txt',
                         'alignment_length', 'perc_similarity',
                         'bit_score', 'by_match', 'per_gene'],
         'orfs': ['initial_summary.txt', 'lengths', 'strands',
                  'by_gene'],
         'ublast_orfs': ['initial_summary.txt',
                         'alignment_length',
                         'perc_similarity',
                         'bit_score',
                         'by_match',
                         'per_gene']}
    outfiles = getOutfiles(D, plot_format)
    return (outfiles)


def makeInitialFigures(n, dpi, by=None):
    figs = []
    if not by:
        by = [4] * n
    for i in np.arange(n):
        figs.append(plt.figure(figsize=(n*by[i], by[i]), dpi=int(dpi)))
    return (figs)


def saveFigures(figs, outfiles, dpi):
    for i, f in enumerate(figs):
        f.savefig(outfiles[i], dpi=int(dpi), bbox_inches='tight')
        plt.close()


def setMax(figure, m, x=True, y=True):
    for ax in figure.get_axes():
        if x:
            ax.set_xlim(min(m['x']), max(m['x'])*1.1)
        if y:
            ax.set_ylim(0, max(m['y'])*1.1)


def roundAll(tab):
    for col in tab.columns:
        if tab[col].dtype == 'float':
            tab[col] = tab[col].round(2)
    return (tab)


def makeSummaryD(table_file,
                 columns, header, columns_to_plot, group_by, maketab=False,
                 keep=""):
    if maketab:
        D = dict()
        F = Fasta.makeFastaDict(table_file)
        for nam in F.keys():
            if keep in nam:
                D[nam] = ORFs.splitNam(nam)
        tab = pd.DataFrame(D).T
    else:
        tab = pd.read_csv(table_file, sep="\t", header=header)
        if columns:
            tab.columns = columns
    if 'length' in columns_to_plot:
        tab['start'] = tab['start'].astype(int)
        tab['end'] = tab['end'].astype(int)
        tab['length'] = tab['end'] - tab['start']
    # calculate the various stats for the summary table
    statD = dict()
    statD['N_Regions_Identified'] = len(tab)
    for col in columns_to_plot:
        statD['Mean_%s' % col] = np.mean(tab[col])
        statD['Max_%s' % col] = np.max(tab[col])
        statD['Min_%s' % col] = np.max(tab[col])

    if group_by:
        grouped = tab.groupby(group_by).size()
        statD['Total_N_%s' % group_by] = len(set(tab[group_by]))
        statD['Mean_by_%s' % group_by] = np.mean(grouped)
        statD['Max_by_%s' % group_by] = np.max(grouped)
        statD['Min_by_%s' % group_by] = np.min(grouped)
        return (statD, tab, grouped)
    return (statD, tab)


def makePlots(plots, cols, colnams, titles, types, table, genes, plotparams,
              outfiles):

    Z = zip(plots, cols, colnams, titles, types)

    for plot, col, colnam, title, typ in Z:
        maxes = {'x': [], 'y': []}
        for i, gene in enumerate(genes):
            sp = plot.add_subplot(1, len(genes), i+1)
            if typ == 'hist':
                h = table[col][table['gene'] == gene]
                sp.hist(h, color=plotparams['%s_colour' % gene])
                sp.set_xlabel(colnam)
                sp.set_ylabel("Frequency")
            elif typ == 'bar':
                table = table.sort_values(col)
                h = table[table['gene'] == gene].groupby(col).size()

                sp.bar(np.arange(len(h)), h,
                       color=plotparams['%s_colour' % gene])
                sp.set_xticks(np.arange(0, len(h)))
                if len(h) > 3:
                    sp.set_xticklabels(h.index, rotation='vertical')
                else:
                    sp.set_xticklabels(h.index)
            sp.set_title(gene)
            maxes['x'].append(sp.get_xlim()[1])
            maxes['x'].append(sp.get_xlim()[0])
            maxes['y'].append(sp.get_ylim()[1])
            maxes['y'].append(sp.get_xlim()[0])
        plot.suptitle(title, y=1.05)
        if plotparams['match_axes'] == "True":
            if typ == 'hist':
                setMax(plot, maxes)
            elif typ == "bar":
                setMax(plot, maxes, x=False)
        plot.tight_layout()
    dpi = int(plotparams['dpi'])
    saveFigures(plots, outfiles, dpi=dpi)


def makeGenePlot(plot, bigtab, genes, title, plotparams, outfile):
    sp = plot.add_subplot(111)
    for i, gene in enumerate(genes):
        sp.bar(i, len(bigtab[bigtab['gene'] == gene]),
               color=plotparams['%s_colour' % gene])
    sp.set_xticks(np.arange(i))
    sp.set_xticklabels(genes)
    sp.set_ylabel("Number of Regions")
    sp.set_title(title)
    dpi = int(plotparams['dpi'])
    saveFigures([plot], [outfile], dpi=dpi)


def summariseExonerateInitial(infiles, outfiles, log, genes, plotparams):
    '''
    Summarise the initial Exonerate run.
    Based on the merged regions in gene_bed_files.dir generated with
    mergeOverlaps.

    '''
    statD = dict()

    # make two empty figures
    dpi = int(plotparams['dpi'])
    lengths, counts, scores, tots, strands = makeInitialFigures(
        5, dpi=dpi)

    log.info("Plotting summary information for initial Exonerate run")
    bigtab = pd.DataFrame()
    biggrouped = pd.DataFrame()
    for i, gene in enumerate(genes):
        # get the files for the right gene
        bed = getInfiles(infiles, gene, 1)
        sD, tab, grouped = makeSummaryD(bed,
                                        HelperFunctions.getBedColumns(),
                                        header=None,
                                        columns_to_plot=['length',
                                                         'bit_score'],
                                        group_by='chrom')
        statD[gene] = sD
        tab['gene'] = gene
        bigtab = bigtab.append(tab)
        grouped = pd.DataFrame(grouped)
        grouped['gene'] = gene
        biggrouped = biggrouped.append(grouped)
    biggrouped.columns = ['count', 'gene']

    makePlots([lengths, scores, strands],
              ['length', 'bit_score', 'strand'],
              ['Length (nt)', 'Exonerate Score', 'Strand'],
              ["Putative ERV Region Lengths - Exonerate Initial",
               "Putative ERV Region Scores - Exonerate Initial",
               "Putative ERV Region Strands - Exonerate Initial"],
              ['hist', 'hist', 'bar'],
              bigtab,
              genes,
              plotparams,
              [outfiles[1], outfiles[2], outfiles[3]])
    makePlots([counts],
              ['count'],
              ['Regions per Sequence'],
              ['Number of ERV Regions Identified per Input Sequence - Exonerate Initial'],
              ['hist'],
              biggrouped,
              genes,
              plotparams,
              [outfiles[4]])

    makeGenePlot(tots, bigtab, genes, "Number of ERV Regions per Gene - Exonerate Initial",
                 plotparams, outfiles[5])
    df = pd.DataFrame(statD)
    df = roundAll(df)
    df.to_csv(outfiles[0], sep="\t")


def summariseUBLAST(infiles, outfiles, log, genes, plotparams):
    dpi = int(plotparams['dpi'])
    tot, lengths, sim, score, target = makeInitialFigures(
        5, dpi=dpi, by=[4, 4, 4, 4, 8])
    log.info("Plotting summary information for UBLAST")
    statD = dict()
    bigtab = pd.DataFrame()
    for i, gene in enumerate(genes):
        f = getInfiles(infiles, gene, 1)
        sD, tab, grouped = makeSummaryD(f,
                                        HelperFunctions.getUBLASTColumns(),
                                        header=None,
                                        columns_to_plot=['alignment_length',
                                                         'percent_identity',
                                                         'bit_score'],
                                        group_by='target')
        tab['target_2'] = ["_".join(x.split("_")[:-3]) for x in tab['target']]
        statD[gene] = sD
        tab['gene'] = gene
        bigtab = bigtab.append(tab)
    makePlots([lengths, sim, score, target],
              ['alignment_length', 'percent_identity', 'bit_score',
               'target_2'],
              ['Length', 'Percent Identity to Target', 'UBLAST Bit Score',
               'Match'],
              ["UBLAST Alignment Lengths of Putative ERV Regions",
               "UBLAST % Identity of Putative ERV Regions",
               "UBLAST Bit Score of Putative ERV Regions",
               "UBLAST Match"],
              ['hist', 'hist', 'hist', 'bar'],
              bigtab,
              genes,
              plotparams,
              outfiles[1:5])
    makeGenePlot(tot, bigtab, genes, "Number of ERV Regions per Gene - UBLAST",
                 plotparams, outfiles[5])
    df = pd.DataFrame(statD)
    df = roundAll(df)
    df.to_csv(outfiles[0], sep="\t")


def FastaLengths(F):
    return (len(x) for x in F.values())


def summariseORFs(infiles, outfiles, log, genes, plotparams):
    dpi = int(plotparams['dpi'])
    tot, lengths, strands = makeInitialFigures(3, dpi=dpi)
    log.info("Plotting summary information for ORFs")

    statD = dict()
    bigtab = pd.DataFrame()
    for i, gene in enumerate(genes):
        f = getInfiles(infiles, gene, 2)

        sD, tab, grouped = makeSummaryD(f,
                                        HelperFunctions.getUBLASTColumns(),
                                        header=None,
                                        columns_to_plot=['length'],
                                        group_by='chrom', maketab=True)
        statD[gene] = sD
        tab['gene'] = gene
        bigtab = bigtab.append(tab)
    makePlots([lengths, strands],
              ['length', 'strand'],
              ['Length', 'strand'],
              ["ORF Length",
               "ORF Strand"],
              ['hist', 'bar'],
              bigtab,
              genes,
              plotparams,
              outfiles[1:3])
    makeGenePlot(tot, bigtab, genes, "Number of ORFs per Gene",
                 plotparams, outfiles[3])
    df = pd.DataFrame(statD)
    df = roundAll(df)
    df.to_csv(outfiles[0], sep="\t")


def summariseGroups(infiles, outfiles, log, genes, plotparams):
    dpi = int(plotparams['dpi'])
    tot, lengths, genus, groups = makeInitialFigures(
        4, dpi=dpi, by=[4, 4, 4, 4, 8])
    log.info("Plotting summary information for Screen section")
    statD = dict()
    bigtab = pd.DataFrame()

    for i, gene in enumerate(genes):
        f = getInfiles(infiles, gene, 0)
        sD, tab, grouped = makeSummaryD(f,
                                        None,
                                        header=0,
                                        columns_to_plot=['length'],
                                        group_by='genus')
        statD[gene] = sD
        tab['gene'] = gene
        bigtab = bigtab.append(tab)
    bigtab = bigtab.sort_values(['gene', 'chrom', 'start'])
    bigtab.to_csv(outfiles[0], sep="\t", index=None)
    makePlots([lengths, genus, groups],
              ['length', 'genus', 'group'],
              ['Length', 'Genus', 'Group'],
              ["Length",
               "Genus",
               "Group"],
              ['hist', 'bar', 'bar'],
              bigtab,
              genes,
              plotparams,
              outfiles[1:4])
    makeGenePlot(tot, bigtab, genes, "Number of ERV Regions per Gene",
                 plotparams, outfiles[4])


def summariseClassify(in_fastas, in_trees, outfiles, genes,
                         plotparams, log):
    log.info("Plotting summary information for groups")
    groups = os.listdir("group_lists.dir")
    rows = []
    for fasta in in_fastas:
        stem = os.path.basename(fasta).split(".")[-2]
        gene, genus = stem.split("_")
        for group in groups:
            if stem in group and group.endswith(".fasta"):
                groupnam = group.split(".")[0]
                F = Fasta.makeFastaDict("group_lists.dir/%s" % group)
                count = len(F)
                rows.append([gene, genus, groupnam, count])
    results = pd.DataFrame(rows, columns=['gene', 'genus', 'group', 'count'])
    results.to_csv(outfiles[0], sep="\t", index=None)
    makeClassifyPlots(results, genes, plotparams, outfiles[1])


def sortkey(x):
    return (x.str.split("_").str.get(-1).astype(int))


def makeClassifyPlots(tab, genes, plotparams, outfile):
    f = plt.figure(figsize=(len(tab) / 10, 12))
    for i, gene in enumerate(genes):

        subtab = tab[tab['gene'] == gene]
        sp = f.add_subplot(3, 1, (i+1))
        x = 0
        splits = []
        counts = []
        ticks = []

        for genus in set(tab['genus']):
            genustab = subtab[subtab['genus'] == genus]
            genustab = genustab.sort_values('group',
                                            key=lambda x: sortkey(x))
            sp.bar(np.arange(x, x+len(genustab)), genustab['count'],
                   color=plotparams['%s_colour' % gene])
            counts += list(genustab['count'])
            ticks += list(genustab['group'])
            x += len(genustab)
            splits.append(x - 0.5)
        sp.vlines(splits, 0, (max(counts) * 1.1))
        sp.set_ylim(0, max(counts) * 1.1)
        sp.set_xlim(0, len(subtab))
        sp.set_xticks(np.arange(len(subtab)))
        sp.set_xticklabels(ticks, rotation='vertical')
        sp.set_xlabel("Group")
        sp.set_ylabel("Number of ORFs")
        sp.set_title("%s %s: %s ORFs" % (gene, genus, sum(subtab['count'])))
    f.tight_layout()
    dpi = int(plotparams['dpi'])
    f.savefig(outfile, dpi=dpi,
              bbox_inches='tight')


def summariseERVRegions(infiles, outfiles, genes, plotparams, log):
    dpi = int(plotparams['dpi'])
    f = plt.figure(figsize=(6, 6), dpi=dpi)
    tab = pd.read_csv(infiles[0], sep="\t")
    tab['typ'] = ["_".join(x.split("_")[:-1]) for x in tab['name']]
    tab = tab.sort_values('typ')
    a = f.add_subplot(111)
    typ_counts = tab.groupby('typ').size()
    a.bar(np.arange(len(typ_counts)), typ_counts,
          color=plotparams['other_colour'])
    a.set_xticklabels(list(typ_counts.index), rotation='vertical')
    a.set_xticks(np.arange(len(typ_counts)))
    a.set_title("Multi-Gene Regions")
    a.set_ylabel("Number of Regions")
    a.set_xlabel("ORFs Present")
    dpi = int(plotparams['dpi'])
    f.savefig(outfiles[1], dpi=dpi, bbox_inches='tight')
