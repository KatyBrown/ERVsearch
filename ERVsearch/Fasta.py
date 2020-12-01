#!/usr/bin/env python3
"""
Functions for reading, generating and processing FASTA files
"""

import os
import pandas as pd
import subprocess
import Errors
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
    # Index the file if required.
    if not os.path.exists("%s.fai" % fasta_in):
        indexFasta(fasta_in, log)

    # Initiate and clear the output file
    if not split:
        out = open(outnam, "w")
        out.close()
    x = 0
    k = 1
    outf = "%s/section1.fasta" % outnam
    for seq in keeplist:
        # get this sequence using samtools faidx
        statement_chrom = ["samtools", "faidx", fasta_in, seq]
        log.info("Processing chromosome %s: %s" % (
            seq, " ".join(statement_chrom)))
        if split:
            if n == 1:
                # if n == 1 write one sequence to each output file and
                # index as you go
                outf = "%s/%s.fasta" % (outnam, seq)
                # put the sequence in its own fasta file
                subprocess.run(statement_chrom, stdout=open(outf, "w"))
                indexFasta(outf, log)
            else:
                if x == n:
                    # if you have enough chromosomes in the current file
                    # index it and make a new file
                    indexFasta(outf, log)

                    outf = "%s/section%i.fasta" % (outnam, k+1)
                    out = open(outf, "w")
                    out.close()
                    k += 1
                    x = 0
                # put the sequence in the ongoing fasta file
                subprocess.run(statement_chrom, stdout=open(outf, "a"))
                x += 1
        else:
            # put the sequence into the one output fasta file
            subprocess.run(statement_chrom, stdout=open(outnam, "a"))

    # index the result if not indexed already
    if not split:
        indexFasta(outnam, log)
    elif n != 1:
        indexFasta(outf, log)


def splitChroms(infile, log, n=1):
    '''
    Split a fasta file into multiple smaller fasta files, with n sequences
    per file.

    If a keep_chroms.txt configuration file is provided, only keep the
    chromosomes listed in this file.
    '''
    indexed = "%s.fai" % os.path.basename(infile)
    allchroms = set([L.strip().split()[0] for L in open(indexed).readlines()])

    keepchroms = readKeepChroms(log)
    if len(keepchroms) == 0:
        keepchroms = allchroms

    # Check all the sequences named in keep_chroms are actually in the input
    # file and error otherwise
    if len(keepchroms & allchroms) != len(keepchroms):
        Errors.raiseError(Errors.MissingChromsError, log=log)

    log.info("Splitting input FASTA file into %i files" % len(keepchroms))
    filterFasta(keepchroms, infile, "host_chromosomes.dir", log, n=n)


def unZip(infile, log):
    '''
    Determines if the file is zipped, gzipped or not zipped and creates
    genome.fa in the working directory as either an unzipped copy or a
    link to the original.
    '''

    if infile.endswith(".gz"):
        # unzip gzipped files
        log.info("Unzipping gzipped input file %s" % infile)
        statement = ["gunzip", "-c", infile]
        log.info("Running statement: %s" % " ".join(statement))
        subprocess.run(statement, stdout=open("genome.fa", "w"))
    elif infile.endswith(".zip"):
        # unzip zipped files
        log.info("Unzipping zipped input file %s" % infile)
        statement = ["unzip", "-p", infile]
        log.info("Running statement: %s" % " ".join(statement))
        subprocess.run(statement, stdout=open("genome.fa", "w"))
    else:
        log.info("Linking to input file %s" % infile)
        statement = ["ln",  "-sf", infile, "genome.fa"]
        subprocess.run(statement)


def indexFasta(infile, log):
    '''
    Indexes a FASTA file using samtools faidx.
    '''
    statement = ["samtools", "faidx", infile]
    log.info("Indexing fasta file: %s" % " ".join(statement))
    s = subprocess.run(statement)
    if s.returncode != 0:
        log.error(s.stderr)
        Errors.raiseError(Errors.FastaIndexError, log=log)


def countFasta(infile, log):
    '''
    Counts the number of chromosomes or sequences in an indexed Fasta file

    '''
    statement = ["wc", "-l", "%s.fai" % infile]
    log.info("Counting chromosomes in fasta file genome.fa: \
             %s" % (" ".join(statement)))
    P = subprocess.run(statement, stdout=subprocess.PIPE)
    nchroms = int(P.stdout.decode().split(" ")[0])
    log.info("%s has %i chromosomes (or contigs / scaffolds" % (infile,
                                                                nchroms))
    return (nchroms)


def readKeepChroms(log):
    '''
    Reads the chromosome names in the keep_chroms.txt file and returns
    them as a set.
    '''
    if os.path.exists("keep_chroms.txt"):
        log.info("""keep_chroms.txt exists so only chromosomes in this\
                    list will be screened""")
        keepchroms = set([L.strip()
                          for L in open("keep_chroms.txt").readlines()])
    else:
        keepchroms = set()
    return (keepchroms)
