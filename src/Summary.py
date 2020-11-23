#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def summariseExonerateInitial(infiles, outfiles, log, genes, plotparams):
    inf = np.array([i[1] for i in infiles])
    inf_stems = np.array([os.path.basename(x) for x in inf])
    statD = dict()
    dpi = int(plotparams['dpi'])
    lengths = plt.figure(figsize=(len(genes)*4, 4), dpi=dpi)
    counts = plt.figure(figsize=(len(genes)*4, 4), dpi=dpi)
    log.info("Plotting summary information for initial Exonerate run")

    for i, gene in enumerate(genes):
        bed = inf[np.char.startswith(inf_stems, gene)][0]

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
    lengths.savefig(outfiles[1], dpi=dpi, bbox_inches='tight')
    counts.savefig(outfiles[2], dpi=dpi, bbox_inches='tight')
