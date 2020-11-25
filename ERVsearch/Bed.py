#!/usr/bin/env python3
"""
Functions for generating and processing bed files
"""

import pandas as pd
import subprocess
pd.set_option('mode.chained_assignment', None)


def combineBeds(beds):
    '''
    Concatenate bed files into a single output and sort
    '''
    rows = []
    for bed in beds:
        with open(bed) as inf:
            for line in inf:
                rows.append(line.strip().split("\t"))

    if len(rows) != 0:
        # sort the combined bed file
        df = pd.DataFrame(rows)
        df[1] = df[1].astype(int)
        df[2] = df[2].astype(int)
        df = df.sort_values([0, 1])
        return (df)


def mergeBed(bed, overlap, log, mergenames=True):
    '''
    Merge overlapping regions in a bed file
    Output the name, strand and score for the highest scoring overlapping
    region.
    '''
    statement = ['bedtools', 'merge',
                 '-s', '-c', '4,5,6', '-d', overlap, '-o', 'collapse']
    log.info("Merging bed file %s : %s" % (bed,
                                           " ".join(statement)))
    P = subprocess.run(statement,
                       stdin=open(bed, "r"),
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.info(P.stderr)
        err = RuntimeError("Bedtools error - see log file")
        log.error(err)
        raise err

    df = pd.DataFrame(
        [x.split("\t") for x in P.stdout.decode().split("\n")[:-1]])
    if mergenames:
        # Pick the highest scoring query for each merged region
        # and output this into the bed file instead of all the names
        rows = []
        for ind in df.index.values:
            row = df.loc[ind]
            scores = [int(x) for x in row[4].split(",")]
            names = row[3].split(",")
            strands = row[5].split(",")
            maxind = scores.index(max(scores))
            maxscore = scores[maxind]
            maxname = names[maxind]
            maxstrand = strands[maxind]
            rows.append([row[0], row[1], row[2],
                         maxname, maxscore, maxstrand])
        merged = pd.DataFrame(rows)
    else:
        merged = df

    return (merged)


def getFasta(infile, outfile, log):
    statement = ['bedtools',
                 'getfasta',
                 '-s', '-fi', "genome.fa",
                 '-bed', infile]
    log.info("Generating fasta file of regions in %s: %s" % (infile,
                                                             statement))

    P = subprocess.run(statement,
                       stdout=open(outfile, "w"),
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr)
        err = RuntimeError("Error converting bed file %s to fasta - \
                           see log file" % infile)
        log.error(err)
        raise err
