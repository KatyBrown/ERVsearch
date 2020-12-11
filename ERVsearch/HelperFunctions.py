#!/usr/bin/env python3
"""
General functions for the pipeline.
"""
import os
import pandas as pd
import shutil
try:
    import ERVsearch.Errors as Errors
except ImportError:
    import Errors

pd.set_option('mode.chained_assignment', None)


def getUBLASTColumns():
    return (['query', 'target', 'percent_identity',
             'alignment_length', 'n_mismatches',
             'n_gap_opens', 'start_pos_query',
             'end_pos_query', 'start_pos_target',
             'end_pos_target', 'evalue', 'bit_score'])


def getBedColumns():
    return (['chrom', 'start', 'end', 'name', 'bit_score', 'strand'])


def revComp(seq):
    '''
    Reverse complements a sequence.
    '''
    rcdict = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N", "Y": "R",
              "K": "M", "R": "Y", "M": "K", "B": "V", "V": "B", "D": "H",
              "H": "D", "W": "W", "S": "S", "-": "-"}
    seq = list(seq)[::-1]
    seq = [rcdict[s] for s in seq]
    seq = "".join(seq)
    return (seq)


def quickCheck(PARAMS, log):
    '''
    Check the environment before starting.
    Check that:
        The input file exists
        The correct path to ERVsearch is provided
        Samtools, Bedtools, FastTree and Mafft are in the PATH
        The correct paths to usearch and Exonerate are provided
    '''
    log.info("Checking ERVsearch setup")
    genome = PARAMS['genome']['file']
    if not os.path.exists(genome):
        if PARAMS['genome']['file'] == "!?":
            Errors.raiseError(Errors.InputFileNotSpecifiedError, log=log)
        else:
            Errors.raiseError(Errors.InputFileNotFoundError, genome,
                              log=log)

    pathD = {'ERVsearch':
             '%s/ERVsearch.py' % (
                 PARAMS['paths']['path_to_ERVsearch']),
             'usearch': PARAMS['paths']['path_to_usearch'],
             'exonerate': PARAMS['paths']['path_to_exonerate']}

    for prog, path in pathD.items():
        if path == "!?":
            Errors.raiseError(Errors.ToolPathNotSpecifiedError, prog, log=log)
        if not os.path.exists(path):
            Errors.raiseError(Errors.ToolNotInPathError, prog, path,
                              log=log)

    for tool in ['samtools', 'bedtools', 'mafft', 'FastTree',
                 'revseq', 'transeq']:
        if not shutil.which(tool):
            Errors.raiseError(Errors.ToolNotInPathError, tool,
                              log=log)
