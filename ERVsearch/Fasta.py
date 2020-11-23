#!/usr/bin/env python3
"""
Functions for reading, generating and processing FASTA files
"""

import os
import pandas as pd
import subprocess
pd.set_option('mode.chained_assignment', None)


def makeFastaDict(multifasta, spliton=None):
    '''
    Turns any fasta file into a dictionary, with sequence names as the keys
    and sequences as the values.
    '''

    querydict = dict()
    x = 0
    seq = []
    nam = ""
    with open(multifasta) as infile:
        for line in infile:
            if line[0] == ">" and x != 0:
                sequ = "".join(seq)
                querydict[nam] = sequ
                seq = []
            if line[0] == ">":
                nam = line.replace(">", "").strip()
                if spliton:
                    nam = nam.split(spliton)[0]
                x += 1
            else:
                seq.append(line.strip())
    sequ = "".join(seq)
    querydict[nam] = sequ
    return querydict


def filterFasta(keeplist, fasta_in, outnam, log, split=True, n=1):
    '''
    Make a FASTA file with only the sequences on keeplist from the input
    fasta file "fasta"
    If split is True and n == 1, make a new fasta file for each sequence
    (with outnam as the directory)

    If split is True and n > 1, put every n sequecnes into a new fasta file
    (with outnam as the directory)

    If split is False put all of the output into the same file
    (with outnam as the file name)
    '''
    if not os.path.exists("%s.fai" % fasta_in):
        statement = ['samtools', 'faidx', fasta_in]
        log.info("Indexing fasta file %s: %s" % (fasta_in,
                                                 " ".join(statement)))
        subprocess.run(statement)
    if not split:
        out = open(outnam, "w")
        out.close()
    x = 0
    k = 1
    outf = "%s/section1.fasta" % outnam
    for seq in keeplist:
        # output a fasta file for each sequence using samtools faidx
        statement_chrom = ["samtools", "faidx", fasta_in, seq]
        log.info("Processing chromosome %s: %s" % (
            seq, " ".join(statement_chrom)))
        if split:
            if n == 1:
                outf = "%s/%s.fasta" % (outnam, seq)
                subprocess.run(statement_chrom, stdout=open(outf, "w"))
                statement_ind = ['samtools', 'faidx', outf]
                log.info("Indexing chromosome %s: %s" % (
                    seq, " ".join(statement_chrom)))
                subprocess.run(statement_ind)
            else:
                if x == n:
                    statement_ind = ["samtools", "faidx", outf]
                    log.info("Indexing %s: %s" % (
                        outf, " ".join(statement_ind)))
                    subprocess.run(statement_ind)

                    outf = "%s/section%i.fasta" % (outnam, k+1)
                    out = open(outf, "w")
                    out.close()
                    k += 1
                    x = 0
                subprocess.run(statement_chrom, stdout=open(outf, "a"))
                x += 1
        else:
            subprocess.run(statement_chrom, stdout=open(outnam, "a"))

    if not split:
        statement = ["samtools", "faidx", outnam]
        log.info("Indexing chromosome %s: %s" % (outnam, " ".join(statement)))
        subprocess.run(statement)
    elif n != 1:
        statement = ["samtools", "faidx", outf]
        log.info("Indexing %s: %s" % (outf, " ".join(statement)))
        subprocess.run(statement)


def splitChroms(infile, log, n=1):
    '''
    For genomes assembled into chromosomes, split the input fasta file into
    one fasta file per chromosome.
    If a keep_chroms.tsv configuration file is provided, only keep the
    chromosomes listed in this file.
    '''
    indexed = "%s.fai" % os.path.basename(infile)
    allchroms = set([L.strip().split()[0] for L in open(indexed).readlines()])

    if os.path.exists("keep_chroms.txt"):
        log.info("""keep_chroms.txt exists so only chromosomes in this\
                    list will be screened""")
        keepchroms = set([L.strip()
                          for L in open("keep_chroms.txt").readlines()])
        if len(keepchroms & allchroms) != len(keepchroms):
            err = RuntimeError("""Not all chromosomes in keepchroms.txt \
                                  are found in the input fasta file, the \
                                  following are missing \n%s""" % (
                                  "\n".join(list(keepchroms - allchroms))))
            log.error(err)
            raise(err)
    else:
        keepchroms = allchroms

    log.info("Splitting input genome into %i chromosomes" % len(keepchroms))
    filterFasta(keepchroms, infile, "host_chromosomes.dir", log, n=n)
