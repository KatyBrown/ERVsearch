#!/usr/bin/env python3
'''
Functions to identify and process open reading frames (ORFs)
'''
import pandas as pd
import numpy as np
import subprocess
import re
import textwrap
try:
    import ERVsearch.Fasta as Fasta
    import ERVsearch.Errors as Errors
except ImportError:
    import Fasta
    import Errors


def runTranseq(fasta, rawout, trans_tab, log):
    '''
    Runs the EMBOSS transeq function to translate a sequence in all six
    frames.
    '''
    # Generate the statement to run transeq
    statement = ['transeq', '-sequence', fasta,
                 '-outseq', rawout,
                 '-frame', '6',
                 '-table', trans_tab,
                 '-auto']
    log.info("Finding ORFs in %s: %s" % (fasta, " ".join(statement)))
    # Run the statement
    P = subprocess.run(statement, stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr)
        Errors.raiseError(Errors.SoftwareError, 'transeq', fasta, log=log)


def splitNam(seqnam):
    '''
    Take a sequence name and return the ID, chromosome, start, end and strand.
    '''
    D = dict()
    groups = re.match(r"([a-z]+[0-9]+)\_(.*)\-([0-9]+)\-([0-9]+)\(([+|-])\)",
                      seqnam)
    D['ID'] = groups[1]
    D['chrom'] = groups[2]
    D['start'] = int(groups[3])
    D['end'] = int(groups[4])
    D['strand'] = groups[5]
    return (D)


def fastaBufferToString(fastabuffer):
    '''
    Converts a fasta file as a buffer read from STDOUT to a string
    with just the sequence.
    '''
    outstr = "".join(fastabuffer.split("\n")[1:])
    return (outstr)


def extractCMD(genome, chrom, start, end, log,
               rc=False, translate=False, trans_table=1):
    '''
    Uses samtools faidx http://www.htslib.org/doc/samtools-faidx.html
    to extract positions start:end from chromosome
    chrom of an indexed genome.
    If rc is True, the sequence is reverse complemented using EMBOSS revcomp
    If translate is True the sequence is translated using the table
    trans_table.
    '''
    # Extract the sequence
    statement = ['samtools', 'faidx', genome,
                 '%s:%i-%i' % (chrom, start, end)]
    log.info("Extracting positions %s to %s from %s:%s : %s" % (
              start, end, genome, chrom, " ".join(statement)))
    P = subprocess.run(statement, stdout=subprocess.PIPE)
    # reverse complement if required
    if rc:
        statement = ['revseq', '-sequence', 'stdin',
                     '-outseq', 'stdout', '-verbose',
                     '0', '-auto']
        log.info("Reverse complementing region %s to %s from %s:%s : %s" % (
                  start, end, genome, chrom, " ".join(statement)))
        P = subprocess.run(statement, stdout=subprocess.PIPE, input=P.stdout,
                           stderr=subprocess.PIPE)

        if P.returncode != 0:
            log.error(P.stderr)
            Errors.raiseError(Errors.SoftwareError,
                              'revseq', "%s %s %s" % (chrom, start, end),
                              log=log)
    out_nt = fastaBufferToString(P.stdout.decode())
    if not translate:
        return (out_nt)
    # translate if required
    statement = ['transeq', '-sequence', 'stdin',
                 '-outseq', 'stdout', '-verbose', '0',
                 '-auto', '-table', str(trans_table)]
    log.info("Translating region %s to %s from %s:%s : %s" % (
              start, end, genome, chrom, " ".join(statement)))
    P = subprocess.run(statement, input=P.stdout,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr)
        Errors.raiseError(Errors.SoftwareError,
                          'transeq', "%s %s %s" % (chrom, start, end),
                          log=log)
    out_aa = fastaBufferToString(P.stdout.decode())
    return (out_nt, out_aa)


def getAllAAs(genome, chrom, cl, start, end, log, trans_table=1):
    '''
    Generates a dictionary with every possible translation of the sequence
    These are extracted directly from the FASTA file.
    '''
    aaD = dict()
    if start == 0:
        # transeq behaves weirdly for position 0
        start += 1
    if end == (cl - 1):
        # and for the last position in the chromosome
        if start > 2:
            start -= 2
        end -= 1
    for i in np.arange(0, 3):
        # extract the region from the chromosome
        nt, aa = extractCMD(genome,
                            chrom,
                            start+i,
                            end,
                            log,
                            trans_table=trans_table,
                            translate=True)
        # store the start position, end position and sequence
        aaD['%s+' % i] = dict()
        aaD['%s+' % i]['aa'] = aa.strip("X")
        aaD['%s+' % i]['start'] = start + i
        aaD['%s+' % i]['end'] = end
        # same for the reverse complement
        nt, aa = extractCMD(genome,
                            chrom,
                            start,
                            end+i,
                            log,
                            rc=True,
                            trans_table=trans_table,
                            translate=True)
        aaD['%s-' % i] = dict()
        aaD['%s-' % i]['aa'] = aa.strip("X")
        aaD['%s-' % i]['start'] = start
        aaD['%s-' % i]['end'] = end + i
    return (aaD)


def filterTranseq(fasta, transeq_raw, nt_out, aa_out, min_length, genome,
                  trans_table, log):
    '''
    Takes the raw output of transeq -frame 6 and finds the longest ORF
    in each region (above a minimum length of min_length).

    ORFs are written to two fasta files - one with the nucleotide sequence
    and one with the amino acid sequence.

    samtools faidx is used to extract the regions from the original genome
    file, EMBOSS revseq is used to reverse complement where needed and
    EMBOSS transeq is used to translate - then these are checked against the
    original ORFs - just to make sure the output is correct.

    ORF names are output as ID_chrom-start-end(strand) with strand relative
    to the original genome input sequences rather than the ORF region
    identified with Exonerate.

    '''
    # get the lengths of the chromosomes - in order to deal correctly with
    # sequence at the ends of the chromosomes
    chrom_df = pd.read_csv("%s.fai" % genome, sep="\t", header=None)
    chrom_lens = dict(zip(chrom_df[0], chrom_df[1]))
    lengths = dict()
    longest = dict()

    # read the raw transeq output into a dictionary
    transeq_raw = Fasta.makeFastaDict(transeq_raw)
    for nam in transeq_raw:
        # split the sequence into ORFs and find the lengths
        Ls = [len(x) for x in transeq_raw[nam].split("*")]
        # if at least one ORF is long enough
        if max(Ls) > min_length:
            log.info("%s has ORFs > %s amino acids: Checking ORFs" % (
                nam, min_length))
            # Get the raw translated sequence generated with transeq
            aa = transeq_raw[nam]

            # The transeq names have _frame as a suffix - cut it off
            nam2 = "_".join(nam.split("_")[:-1])
            lengths.setdefault(nam2, 0)
            longest.setdefault(nam2, dict())

            # extract the ID, start and end positions and strand from the
            # sequence name
            namD = splitNam(nam2)

            # get every possible translation by manually RCing and clipping
            # the sequences
            aaD = getAllAAs(genome, namD['chrom'], chrom_lens[namD['chrom']],
                            namD['start'], namD['end'], log,
                            trans_table=trans_table)

            # iterate through the translation frames
            for f, new_aa in aaD.items():
                # this should be the correct translation
                if aa in new_aa['aa']:
                    if "+" in f:
                        # adjust so that the ends are the same as the
                        # original translation
                        mod_s = new_aa['aa'].find(aa) * 3
                        mod_e = new_aa['aa'][::-1].find(aa[::-1]) * 3
                        s = new_aa['start'] + mod_s
                        e = new_aa['end'] - mod_e
                        rc = False
                    else:
                        mod_s = new_aa['aa'].find(aa) * 3
                        mod_e = new_aa['aa'][::-1].find(aa[::-1]) * 3
                        s = new_aa['start'] + mod_e - 1
                        e = new_aa['end'] - mod_e
                        rc = True
                    # split the amino acid sequence into ORFs
                    orfs = aa.split("*")
                    for orf in orfs:
                        if len(orf) > 50:
                            opos = aa.find(orf)
                            o_nt_pos = (opos * 3)
                            if "+" in f:
                                orf_s = s + o_nt_pos
                                orf_e = s + o_nt_pos + (len(orf) * 3) - 1
                            else:
                                orf_s = e - o_nt_pos - (len(orf) * 3)
                                orf_e = e - o_nt_pos

                            # extract the ORF directly from the FASTA file
                            nt, aa2 = extractCMD(genome,
                                                 namD['chrom'],
                                                 orf_s,
                                                 orf_e,
                                                 log,
                                                 rc=rc,
                                                 trans_table=trans_table,
                                                 translate=True)
                            aa2 = aa2.strip("X")
                            orf = orf.strip("X")
                            # Adjust again - sometimes there are stop codons
                            if aa2.startswith("*"):
                                orf_s -= 3
                                orf_e -= 3
                            if aa2.endswith("*"):
                                orf_s += 3
                                orf_e += 3
                            if orf_e > chrom_lens[namD['chrom']]:
                                orf_e = chrom_lens[namD['chrom']]
                            nt, aa2 = extractCMD(genome,
                                                 namD['chrom'],
                                                 orf_s,
                                                 orf_e,
                                                 log,
                                                 rc=rc,
                                                 trans_table=trans_table,
                                                 translate=True)
                            # Remove Xs from the ends of the sequences
                            aa2 = aa2.strip("X")
                            orf = orf.strip("X")
                            if "+" in f:
                                frame = "+"
                            else:
                                frame = "-"
                            # Create a new ID
                            ID = "%s_%s-%s-%s(%s)" % (namD['ID'],
                                                      namD['chrom'],
                                                      orf_s,
                                                      orf_e,
                                                      frame)
                            if len(aa2) > lengths[nam2]:
                                lengths[nam2] = len(aa2)
                                longest[nam2]['aa'] = aa2
                                longest[nam2]['nt'] = nt
                                longest[nam2]['ID'] = ID
    nt_out = open(nt_out, "w")
    aa_out = open(aa_out, "w")
    # output the longest ORF in each ERV region to file
    for nam in longest:
        if 'ID' in longest[nam]:
            nt = "\n".join(textwrap.wrap(longest[nam]['nt'], 70))
            aa = "\n".join(textwrap.wrap(longest[nam]['aa'], 70))
            nt_out.write(">%s\n%s\n" % (longest[nam]['ID'], nt))
            aa_out.write(">%s\n%s\n" % (longest[nam]['ID'], aa))
    nt_out.close()
    aa_out.close()
