#!/usr/bin/env python3
"""
Functions for generating summary tables and plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def getInfiles(allinfiles, gene, ind):
    '''
    Get the input files corresponding to a particular index and gene
    '''
    inf = np.array([i[ind] for i in allinfiles])
    inf_stems = np.array([os.path.basename(x) for x in inf])
    infile = inf[np.char.startswith(inf_stems, gene)][0]
    return (infile)


def makeInitialFigures(n, dpi):
    figs = []
    for i in np.arange(n):
        figs.append(plt.figure(figsize=(n*4, 4), dpi=int(dpi)))
    return (figs)


def saveFigures(figs, outfiles, dpi):
    for i, f in enumerate(figs):
        f.savefig(outfiles[i], dpi=int(dpi), bbox_inches='tight')


def summariseExonerateInitial(infiles, outfiles, log, genes, plotparams):
    '''
    Summarise the initial Exonerate run.
    '''
    statD = dict()

    lengths, counts = makeInitialFigures(2, dpi=plotparams['dpi'])

    log.info("Plotting summary information for initial Exonerate run")

    for i, gene in enumerate(genes):
        bed = getInfiles(infiles, gene, 1)

        bed_tab = pd.read_csv(bed, sep="\t", header=None)
        bed_tab.columns = ['chromosome', 'start', 'end', 'best_match',
                           'exonerate_score', 'strand']

        bed_tab['length'] = bed_tab['end'] - bed_tab['start']

        by_chrom = bed_tab.groupby('chromosome').size()
        L = bed_tab['length']
        statD[gene] = {"N_Regions_Identified": len(bed_tab),
                       "Mean_Length_Regions": np.mean(L),
                       'Max_Length_Regions': np.max(L),
                       'Min_Length_Regions': np.min(L),
                       "N_Input_Sequences_With_Regions": len(set(
                           bed_tab['chromosome'])),
                       'Mean_Regions_Per_Input_Sequence': np.median(by_chrom),
                       'Maximum_Regions_Per_Input_Sequence': np.max(by_chrom),
                       'Minimum_Regions_Per_Input_Sequence': np.min(by_chrom)}
        ls = lengths.add_subplot(1, len(genes), i+1)
        ls.hist(bed_tab['length'], color=plotparams['%s_colour' % gene])
        ls.set_xlabel("Length (nt)")
        ls.set_ylabel("Frequency")
        ls.set_title(gene)
        c = counts.add_subplot(1, len(genes), i+1)
        c.hist(by_chrom, color=plotparams['%s_colour' % gene])
        c.set_xlabel("Count per Sequence")
        c.set_ylabel("Frequency")
        c.set_title(gene)
    lengths.suptitle("Putative ERV Region Lengths - Exonerate Initial",
                     y=1.1)
    counts.suptitle(
        "Number of Putative ERV Regions per Sequence - Exonerate Initial",
        y=1.1)

    lengths.tight_layout()
    counts.tight_layout()
    df = pd.DataFrame(statD)
    df.to_csv(outfiles[0], sep="\t")
    saveFigures([lengths, counts], outfiles[1:], plotparams['dpi'])


def summariseUBLAST(infiles, outfiles, log, genes, plotparams):
    lengths, sim, score = makeInitialFigures(3, dpi=plotparams['dpi'])
    log.info("Plotting summary information for UBLAST")
    for gene in genes:
        f = getInfiles(infiles, gene, 1)
    
        tab = pd.read_csv(bed, sep="\t", header=None)
        tab.columns = ['chromosome', 'start', 'end', 'best_match',
                           'exonerate_score', 'strand']
    
        bed_tab['length'] = bed_tab['end'] - bed_tab['start']
    
        by_chrom = bed_tab.groupby('chromosome').size()
        L = bed_tab['length']
        statD[gene] = {"N_Regions_Identified": len(bed_tab),
                       "Mean_Length_Regions": np.mean(L),
                       'Max_Length_Regions': np.max(L),
                       'Min_Length_Regions': np.min(L),
                       "N_Input_Sequences_With_Regions": len(set(
                           bed_tab['chromosome'])),
                       'Mean_Regions_Per_Input_Sequence': np.median(by_chrom),
                       'Maximum_Regions_Per_Input_Sequence': np.max(by_chrom),
                       'Minimum_Regions_Per_Input_Sequence': np.min(by_chrom)}
        ls = lengths.add_subplot(1, len(genes), i+1)
        ls.hist(bed_tab['length'], color=plotparams['%s_colour' % gene])
        ls.set_xlabel("Length (nt)")
        ls.set_ylabel("Frequency")
        ls.set_title(gene)
        c = counts.add_subplot(1, len(genes), i+1)
        c.hist(by_chrom, color=plotparams['%s_colour' % gene])
        c.set_xlabel("Count per Sequence")
        c.set_ylabel("Frequency")
        c.set_title(gene)    