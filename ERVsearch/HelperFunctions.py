#!/usr/bin/env python3
"""
General functions for the pipeline.
"""
import os
import pandas as pd
import shutil

pd.set_option('mode.chained_assignment', None)


def getUBLASTColumns():
    return (['query', 'target', 'percent_identity',
             'alignment_length', 'n_mismatches',
             'n_gap_opens', 'start_pos_query',
             'end_pos_query', 'start_pos_target',
             'end_pos_target', 'evalue', 'bit_score'])


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
        samtools, bedtools, FastTree and mafft are in the PATH
        The correct paths to usearch and exonerate are provided
    '''
    if not os.path.exists(PARAMS['genome']['file']):
        if PARAMS['genome']['file'] == "!?":
            err = RuntimeError(
                "Input file needs to be specified in the pipeline.ini file")
        else:
            err = FileNotFoundError("Input file %s not found" % (
                PARAMS['genome']['file']))
        log.error(err)
        raise (err)

    pathD = {'ERVsearch':
             '%s/src/pipeline_ERVs.py' % (
                 PARAMS['database']['path_to_ERVsearch']),
             'usearch': PARAMS['paths']['path_to_usearch'],
             'exonerate': PARAMS['paths']['path_to_exonerate']}

    for prog, path in pathD.items():
        if not os.path.exists(path):
            err = RuntimeError("%s is not at the location %s specified in \
                               your pipeline.ini" % (prog, path))
            log.error(err)
            raise (err)

    for tool in ['samtools', 'bedtools', 'mafft', 'FastTree',
                 'revseq', 'transeq']:
        if not shutil.which(tool):
            err = RuntimeError("%s is not in your $PATH" % tool)
            log.error(err)
            raise(err)
