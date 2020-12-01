#!/usr/bin/env python3
"""
Functions for generating and processing bed files
"""

import pandas as pd
import subprocess
pd.set_option('mode.chained_assignment', None)


def combineBeds(beds):
    '''
    Concatenate bed files into a single output and sort by chromosome then
    start position
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


def mergeBeds(bed, overlap, log, mergenames=True):
    '''
    Merge overlapping regions in a bed file

    If mergenames is True, then for merged regions columns 4-6 will contain
    the name, score and strand of the region which was merged in which
    has the highest score
    (based on whichever score is in column 5 of the input)

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
            # get all the scores
            scores = [int(x) for x in row[4].split(",")]
            names = row[3].split(",")
            strands = row[5].split(",")

            # find the index of the highest scoring merged region
            maxind = scores.index(max(scores))
            # find the value of the highest scoring region
            maxscore = scores[maxind]
            # name of the highest scoring region
            maxname = names[maxind]
            # strand of the highest scoring region
            maxstrand = strands[maxind]
            # put it into the table
            rows.append([row[0], row[1], row[2],
                         maxname, maxscore, maxstrand])
        merged = pd.DataFrame(rows)
    else:
        merged = df

    return (merged)


def getFasta(infile, outfile, log):
    '''
    Extracts the regions specified in the genome from a FASTA file
    using bedtools getfasta
    https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html

    The genome has to be saved as genome.fa in the working directory (the
    pipeline sets this up).
    '''
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
