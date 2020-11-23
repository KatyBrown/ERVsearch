#!/usr/bin/env python3
"""
Functions to run Exonerate and process the output.
"""

import os
import pandas as pd
import numpy as np
import subprocess
import copy


def runExonerate(fasta, chrom, outf, log, exonerate):
    '''
    Runs Exonerate protein2dna.
    Query - FASTA file containing retrovirus amino acid sequences.

    Target - chromosomes in the host_chromosomes directory.
    Settings have been optimised for time and ERV detection.
    '''
    if not os.path.exists(exonerate):
        err = FileNotFoundError("The path %s to the exonerate executable \
                                is not valid" % exonerate)
        log.error(err)
        raise err
    statement = [exonerate, '--model', 'protein2dna',
                 '--showalignment', 'F', '--seedrepeat', '1',
                 '--showvulgar', 'T', '--query', fasta,
                 '--target', chrom]
    log.info("Running Exonerate on %s vs %s: %s" % (
        chrom, fasta, " ".join(statement)))
    P = subprocess.run(statement, universal_newlines=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stdout)
        err = RuntimeError("Exonerate error - see log file")
        log.error(err)
        raise err
    else:
        out = open(outf, "w")
        out.write(P.stdout)
        out.close()


def filterExonerate(infile, outfiles, min_hit_length, log):
    '''
    Converts the results of the exonerate screen into a pandas dataframe
    and saves as XXX_unfiltered.tsv.
    Filters hits shorter than the min_hit_length specified in pipeline.ini
    and hits containing introns and saves as XXX_filtered.tsv
    Converts to bed format - XXX.bed
    '''
    # Clean column names.  Details are from the final column of the
    # Exonerate showvulgar output, as listed here
    # http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual

    cnames = ["query_id", "query_start", "query_end", "query_strand",
              "target_id", "target_start", "target_end",
              "target_strand", "score", "details"]
    log.info("Converting Exonerate output %s to dataframe" % outfiles[0])
    rows = []
    with open(infile) as inf:
        for line in inf:
            if line.startswith('vulgar'):
                line = line.strip().split(" ")
                details = "|".join(line[10:])
                rows.append(line[1:10] + [details])
    log.info("%i rows identified in %s" % (len(rows), infile))
    exonerate = pd.DataFrame(rows, columns=cnames)
    for cname in ['query_start', 'query_end', 'target_start', 'target_end',
                  'score']:
        exonerate[cname] = exonerate[cname].astype(int)
    exonerate['length'] = abs(exonerate['target_end'] -
                              exonerate['target_start'])

    # swap the target start and target end co-ordinates for minus strand
    # otherwise bedtools getfasta doesn't work
    minus = copy.copy(exonerate[exonerate['target_strand'] == "-"])
    exonerate['target_start'][exonerate['target_strand'] == "-"] = minus[
        'target_end']
    exonerate['target_end'][exonerate['target_strand'] == "-"] = minus[
        'target_start']
    exonerate = exonerate.sort_values('target_start')
    exonerate.to_csv(outfiles[0], sep="\t", index=None)

    log.info("Filtering Exonerate output %s" % outfiles[0])
    # Filter hits which are too small or contain introns
    exonerate = exonerate[((exonerate['length'] > min_hit_length) &
                           (~exonerate['details'].str.contains("I")))]
    log.info("%i rows filtered out of %s" % (len(rows) - len(exonerate),
                                             outfiles[1]))
    exonerate.to_csv(outfiles[1], sep="\t", index=False)

    log.info("Converting Exonerate output %s to bed format" % outfiles[1])
    bedcols = exonerate[['target_id', 'target_start', 'target_end', 'query_id',
                         'score', 'target_strand']]
    bedcols.to_csv(outfiles[2],
                   sep="\t", header=False, index=False)


def classifyWithExonerate(reference_ervs, fasta, raw_out,
                          exonerate, exonerate_minscore,
                          log):
    '''
    Runs the exonerate ungapped algorithm with each ERV region
    in the fasta file against the
    all_ERVS.fasta fasta file, to detect which known
    retrovirus is most similar to each newly identified ERV region.
    '''

    # we can use a more simple output format this time because we're not
    # interested in introns

    # generate the statement to run exonerate
    statement = [exonerate,
                 '--query', fasta,
                 '--target', reference_ervs,
                 '--showalignment', 'F',
                 '--showvulgar', 'F',
                 '--ryo', "%qi\\t%ti\\t%s\\n",
                 '--score', exonerate_minscore]
    log.info("Running a second exonerate pass for classification \
              on %s: %s" % (fasta, " ".join(statement)))

    P = subprocess.run(statement, universal_newlines=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)

    if P.returncode != 0:
        log.error(P.stdout)
        err = RuntimeError("Exonerate error - see log file")
        log.error(err)
        raise err
    else:
        out = open(raw_out, "w")
        out.write(P.stdout)
        out.close()


def findBestExonerate(res, gene):
    '''
    Takes the output of running Exonerate on all putative ERVs
    vs a database of known ERVs, filters to keep only those matching the
    correct gene then picks the highest scoring result for each putative
    ERV.
    '''
    res['gene'] = res['match'].str.split("_").str.get(-1)
    # check the hits are for the right gene
    res = res[res['gene'] == gene]
    res = res.drop('gene', 1)
    # Sort by score
    res = res.sort_values('score', ascending=False)
    # Take the highest scoring result for each ID
    res = res.groupby('id').first()
    res['id'] = res.index.values
    res.index = np.arange(len(res))
    res['genus'] = res['match'].str.split("_").str.get(-2)
    return (res)
