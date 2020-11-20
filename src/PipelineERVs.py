# This file provides functions used by the pipeline.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
import ete3 as ete
import subprocess
import copy
import shutil
import re
import textwrap

pd.set_option('mode.chained_assignment', None)

cols = ['crimson', 'mediumspringgreen', 'deepskyblue',
        'goldenrod', 'deeppink', 'mediumpurple', 'orangered']
genera = ['gamma', 'beta', 'spuma', 'epsilon', 'alpha', 'lenti', 'delta']
coldict = dict(zip(genera, cols))


def zapFile(filename):
    '''
    Removes the contents of a file but maintains the time stamp

    Used to remove input files after output has been generated but not
    disrupt the pipeline.
    '''
    # store the time stamp from the file
    original = os.stat(filename)
    # remove file contents
    f = open(filename, "w")
    f.close()
    # change time stamp to original
    os.utime(filename, (original.st_atime, original.st_mtime))


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


def mergeBed(bed, overlap, log):
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

    # Pick the highest scoring query for each merged region
    # and output this into the bed file instead of all
    # the names
    df = pd.DataFrame([
        x.split("\t") for x in P.stdout.decode().split("\n")[:-1]])
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
    return (merged)


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


def simpleFasta(L, seqs, outfile):
    '''
    Generates a fasta file containing sequences in list L, a subset of a
    larger fasta file seqs.
    '''

    fastaD = makeFastaDict(seqs)
    outf = open(outfile, "w")
    for line in L:
        seq = fastaD[line]
        outf.write(">%s\n%s\n" % (line, seq))
    outf.close()


def keyfunc(s):
    '''
    Helper function for natural sort of chromosomes by number rather than
    name
    '''

    return [int(''.join(g)) if k else ''.join(g)
            for k, g in groupby(s, str.isdigit)]


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


def runTranseq(fasta, rawout, trans_tab, log):

    if not shutil.which('transeq'):
        err = FileNotFoundError("EMBOSS transeq is not in your PATH")
        log.error(err)
        raise err
    statement = ['transeq', '-sequence', fasta,
                 '-outseq', rawout,
                 '-frame', '6',
                 '-table', trans_tab,
                 '-auto']
    log.info("Finding ORFs in %s: %s" % (fasta, " ".join(statement)))
    P = subprocess.run(statement, stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr)
        err = RuntimeError(
            "Error running sixpack on %s - see log file" % fasta)
        log.error(err)
        raise err


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
    Uses samtools faidx to extract positions start:end from chromosome
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
        P = subprocess.run(statement, stdout=subprocess.PIPE, input=P.stdout)
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
                       stdout=subprocess.PIPE)

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

        nt, aa = extractCMD(genome,
                            chrom,
                            start+i,
                            end,
                            log,
                            trans_table=trans_table,
                            translate=True)
        aaD['%s+' % i] = dict()
        aaD['%s+' % i]['aa'] = aa.strip("X")
        aaD['%s+' % i]['start'] = start + i
        aaD['%s+' % i]['end'] = end
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
    chrom_df = pd.read_csv("%s.fai" % genome, sep="\t", header=None)
    # get the lengths of the chromosomes
    chrom_lens = dict(zip(chrom_df[0], chrom_df[1]))
    lengths = dict()
    longest = dict()
    transeq_raw = makeFastaDict(transeq_raw)
    for nam in transeq_raw:
        Ls = [len(x) for x in transeq_raw[nam].split("*")]
        if max(Ls) > min_length:
            log.info("%s has ORFs > %s amino acids: Checking ORFs" % (
                nam, min_length))
            # Get the raw translated sequence generated with transeq
            aa = transeq_raw[nam]
            nam2 = "_".join(nam.split("_")[:-1])
            lengths.setdefault(nam2, 0)
            longest.setdefault(nam2, dict())

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
                            aa2 = aa2.strip("X")
                            orf = orf.strip("X")
                            if "+" in f:
                                frame = "+"
                            else:
                                frame = "-"
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
    for nam in longest:
        if 'ID' in longest[nam]:
            nt = "\n".join(textwrap.wrap(longest[nam]['nt'], 70))
            aa = "\n".join(textwrap.wrap(longest[nam]['aa'], 70))
            nt_out.write(">%s\n%s\n" % (longest[nam]['ID'], nt))
            aa_out.write(">%s\n%s\n" % (longest[nam]['ID'], aa))
    nt_out.close()
    aa_out.close()


def group(infile, convert, outfile):
    '''
    The convert file approximately groups known retroviruses by phylogenetic
    group.  This function uses the output of the findBest function in
    the pipeline file to assign newly identified ERVs to these
    approximate groups.
    '''

    try:
        df = pd.read_csv(infile, sep="\t", index_col=0)
    except:
        df = []
    if len(df) != 0:
        # Expands the sequence names to get the start and end positions of
        # the ERV like regions as the ORF step does not output these.  Merges
        # these back into the input and converts them to integers.
        details = df['id'].str.split(":", expand=True)
        details.columns = ['chr', 'se']
        details2 = details['se'].str.split("-", expand=True)
        details2.columns = ['start', 'end']
        details = details.drop('se', 1)
        details = details.merge(details2, left_index=True, right_index=True)
        details['start'] = details['start'].astype(int)
        details['end'] = details['end'].astype(int)
        details['length'] = details['end'] - details['start']

        # Find the group in convert of the info in the "match"
        # column in the dataframe and  add this to the "group" column
        df['group'] = df['match'].apply(
            lambda x: x if x not in convert['id']
            else convert['match'][convert['id'] == x].values[0])
        res = df.merge(details, left_index=True, right_index=True)
        res.to_csv(outfile, sep="\t")
    else:
        # If there are no hits for this gene, generate a blank outfile.
        os.system("touch %s" % outfile)


def makeNovelDict(grouptable, group, novelfasta):
    '''
    Used by makePhyloFasta
    Builds a dictionary of the names and sequences of the newly identified
    ERV regions.
    If there are more than 40 sequences in a group, a random sample of 20
    sequences is extracted (as very large trees tend to be inaccurate)
    '''

    noveldict = dict()
    groupseqs = grouptable[grouptable['group'] == group]

    idlist = groupseqs['id']
    idlist = np.unique(idlist)
    for id in idlist:
        noveldict[id] = novelfasta[id]
    return noveldict


def makeRefDict(phylopath, gene, genus, group):

    '''
    Used by makePhyloFasta
    Takes a set of known retrovirus sequences corresponding to the group
    of interest and builds a dictionary where keys are names and values
    are sequences.
    If there are few or no known sequences in the group, a standard set
    of reference sequences are used.
    '''

    reffasta = makeFastaDict("%s/%s_%s.fasta" % (phylopath, gene, genus))
    try:
        groupfasta = makeFastaDict("%s/%s.fasta" % (phylopath, group))
    except:
        groupfasta = makeFastaDict("%s/%s_%s_strays.fasta"
                                   % (phylopath, gene, genus))
        groupfasta.update(reffasta)

    if group not in groupfasta.keys():
        groupfasta = reffasta

    # If there are less than 8 sequences in the group, add the standard
    # reference sequences (to avoid very small trees)
    if len(groupfasta) < 8:
        groupfasta.update(reffasta)

    return groupfasta


def getOutgroup(phylopath, gene, genus):

    '''
    Uses the provided reference file to identify an appropriate outgroup
    for a particular gene and genus of ERVs.
    '''

    outgdict = [tuple(line.strip().split("\t"))
                for line in open(
                        "%s/outgroups.tsv" % phylopath).readlines()]
    outgdict = dict(outgdict)
    genegenus = "%s_%s" % (genus, gene)
    outg = outgdict[genegenus]
    return outg


def makePhyloFasta(grouptable, allgroups, fasta, phylopath, mafftpath,
                   fasttreepath, seqdir):

    '''
    Builds a phylogenetic tree representing each small group of ERVs
    Groups are identified in the makeGroups stage.
    Sequences are compiled into fasta files, aligned using MAFFT and
    trees are built using FastTree. An image of each tree is also generated.
    '''

    novelfasta = makeFastaDict(fasta)
    allfasta = makeFastaDict("%s/all_ERVS.fasta" % seqdir)
    for group in allgroups:
        gene = group.split("_")[-1]
        genus = group.split("_")[-2]

        # Build dictionaries of sequences and sequence names
        iddict = makeNovelDict(grouptable, group, novelfasta)
        refdict = makeRefDict(phylopath, gene, genus, group)

        outg = getOutgroup(phylopath, gene, genus)
        refdict[outg] = allfasta[outg]

        #  write the fasta files
        outf = ("fastas/%s.fasta" % (group))
        out = open(outf, "w")

        for key in refdict:
            seq = refdict[key]
            out.write(">%s\n%s\n" % (key, seq))

        for id in iddict:
            seq = iddict[id]
            out.write(">%s\n%s\n" % (id.replace(":", "_"), seq))

        out.close()

        # Align Fasta Files
        ali = outf.replace(".fasta", "_ali.fasta")
        tree = ali.replace("_ali.fasta", ".tre")
        tree = tree.replace("fastas", "trees")
        os.system("fftns --quiet %s > %s"
                  % (outf, outf.replace(".fasta", "_ali.fasta")))

        # Build Trees
        os.system("%s -nt -gtr -quiet %s > %s" % (fasttreepath,
                                                  ali, tree))

        # Build Tree Images
        treepng = tree.replace(".tre", ".png")
        T = ete.Tree(tree)
        for item in T.traverse():
            if item.name[0:3] == "chr":
                col = 'darkgreen'
            elif item.name == outg:
                col = 'blue'
            else:
                col = 'black'

            TF = ete.TextFace(item.name)
            TF.fgcolor = col
            item.add_face(TF, column=0)
        T.set_outgroup(outg)
        TS = ete.TreeStyle()
        TS.show_leaf_name = False
        T.render(treepng, tree_style=TS)


def remove_subsets(set1, list_of_sets):

    '''
    Used by monophyleticERVgroups
    Takes a specific set of strings and a list of sets of strings.
    Sets in the list are removed if they are subsets of the specific set.
    '''

    i = 0
    L = []
    for aset in list_of_sets:
        if not aset < set1:
            L.append(i)
        i += 1
    L = np.array(L)
    S = np.array(list_of_sets)
    return S[L]


def monophyleticERVgroups(tree):

    '''
    Based on a phylogenetic tree, finds monophyletic groups containing only
    novel sequences (identified because names start with gag, pol or env)
    '''
    p = ete.Tree(tree)
    sets = []
    done = set()
    genes = ['gag', 'pol', 'env']
    # Runs through all the subtrees of the specified tree and determines
    # if they contain any non-novel groups.
    for x in p.traverse():
        leafnames = x.get_leaf_names()
        count = 0
        for leaf in leafnames:
            if leaf[0:3] in genes and leaf not in done:
                count += 1
        if count == len(leafnames):
            done = done | set(leafnames)
            sets.append(set(leafnames))

    # sort the groups by desending size (by number of elements)
    return sets



def makeRepFastas(infiles, pp,
                  mafftpath, fasttreepath, seqdir, outfile):
    '''
    Using the Fasta files from makePhyloFastas, takes a single representative
    for monophyletic groups and builds a summary tree for each genus and gene.
    Generates a phylogeny and an image of each phylogeny with node labels sized
    according to the size of the group.
    '''
    genus = outfile.split("/")[-1].split("_")[-2]
    gene = outfile.split("_")[-1][0:3]

    allfasta = makeFastaDict("%s/all_ERVS.fasta" % seqdir)
    repcounts = dict()
    repseqs = dict()
    repgroups = dict()

    # Choose a representative for each monophyletic cluster of novel ERV
    # regions in the input trees
    for tree in infiles:
        fasta = tree.replace("tree", "fasta")
        fasta = fasta.replace(".tre", ".fasta")
        fasta = makeFastaDict(fasta, spliton="(")

        stem = tree.split("/")[-1].split(".")[0]
        groups = monophyleticERVgroups(tree)

        # Record the representative member of each group, its sequence and
        # the corrected group size.  If the input trees contain all
        # sequences, rather than subsets of sequences, make an outfile
        # with an ID number containing all the sequence names in the group

        j = 0
        k = 1
        for group in groups:
            L = list(group)
            rep = L[0]
            repcounts[rep] = len(groups[j])
            repseqs[rep] = fasta[rep]

            groupid = "%s_%s_%i" % (genus, gene, k)
            gout = open("group_lists/%s.txt" % groupid, "w")
            gfasta = open("group_fastas/%s.fasta" % groupid, "w")
            for nam in group:
                gfasta.write(">%s\n%s\n" % (nam, fasta[nam]))
                if nam == rep:
                    gout.write("%s**\n" % nam)
                else:
                    gout.write("%s\n" % nam)
            gout.close()
            gfasta.close()
            repgroups[rep] = groupid
            j += 1
            k += 1

    # Build a fasta file with the group representatives and the representative
    # reference trees

    others = makeFastaDict("%s/%s_%s.fasta" % (pp, gene, genus))
    outg = getOutgroup(pp, gene, genus)
    others[outg] = allfasta[outg]

    fastadir = outfile.replace("summary_trees", "summary_fastas")

    fastaout = fastadir.replace(".tre", ".fasta")
    aliout = fastadir.replace(".tre", "_ali.fasta")

    out = open(fastaout, "w")
    for other in others:
        out.write(">%s\n%s\n" % (other, others[other]))

    for novel in repseqs:
        out.write(">%s\n%s\n" % (novel, repseqs[novel]))
    out.close()

    # Align the Fasta file
    os.system("fftns --quiet %s > %s"
              % (fastaout,
                 fastaout.replace(".fasta", "_ali.fasta")))

    # Build a tree of the fasta file
    os.system("FastTree -nt -gtr -quiet %s > %s"
              % (aliout, outfile))

    # Build an image of the fasta file
    treepng = outfile.replace(".tre", ".png")
    T = ete.Tree(outfile)
    m = max(repcounts.values())
    scale = 40 / float(m)

    for item in T.traverse():
        if item.name[0:3] in ['gag', 'pol', 'env']:
            group = repgroups[item.name]
            count = repcounts[item.name]
            countn = round(scale * count, 0)
            C = ete.CircleFace(radius=countn, color=coldict[genus])
            item.add_face(C, column=0)
            R = ete.TextFace("  %s" % count)
            R.fgcolor = 'black'
            R.fsize = 16
            item.add_face(R, column=1)
            TF = ete.TextFace("  %s" % group)
            TF.fgcolor = 'black'
            item.add_face(TF, column=2)
        elif item.name == outg:
            nam = " ".join(item.name.split("_")[0:5])
            nam = nam.replace("Human_ERV_", "HERV_")
            TF = ete.TextFace(nam)
            TF.fgcolor = 'blue'
            TF.fsize = 12
            item.add_face(TF, column=1)
        else:
            nam = " ".join(item.name.split("_")[0:5])
            nam = nam.replace("Human_ERV_", "HERV_")
            TF = ete.TextFace(nam)
            TF.fsize = 12
            TF.fgcolor = 'black'
            item.add_face(TF, column=1)

    T.set_outgroup(outg)
    TS = ete.TreeStyle()
    titl = ete.TextFace("%s, %s,  %i seqs"
                        % (gene, genus, sum(repcounts.values())),
                        fsize=20)
    titl.margin_bottom = 20
    TS.title.add_face(titl, column=1)
    TS.show_leaf_name = False
    TS.show_branch_support = True
    T.render(treepng, tree_style=TS, dpi=600)


def summary(gag, pol, env, outfile, plot_dir):
    '''
    Using matplotlib, builds a series of summary plots about the ERVs
    detected for each gene.
    '''
    genes = ['Gag', 'Pol', 'Env']
    out = open(outfile, "w")
    generaL = pd.Series(genera)
    cols = [coldict[genus] for genus in generaL]
    i = 0
    for parsed_output in (gag, pol, env):
        # Write basic statistics to a file - summary.tsv
        gene = genes[i]
        if len(parsed_output) == 0:
            out.write("No results for %s\n" % gene)
        else:

            out.write("Number of Hits %s = %i\n" % (gene,
                                                    len(parsed_output)))
            M = np.max(parsed_output['orf_len'])
            out.write("Longest ORF %s = %i\n" % (gene, M))

            # Plot a histogram of the ORF length
            f = plt.figure()
            a = f.add_subplot('111')
            a.hist(parsed_output['orf_len'], bins=20, color="crimson")
            a.xaxis.set_label_text('ORF length')
            a.yaxis.set_label_text('Frequency')
            f.suptitle('ORF Length %s' % gene)
            f.savefig("%s/orf_lengths_%s.png" % (plot_dir, gene),
                      bbox_inches='tight')

            # Plot a bar chart of the number of ERVs per chromosome
            counts = np.unique(parsed_output['chr'], return_counts=True)
            i2 = sorted(counts[0], key=keyfunc)
            cframe = pd.DataFrame(counts[1], counts[0],
                                  columns=['counts']).reindex(index=i2)

            fig = plt.figure()
            sp = fig.add_subplot('111')
            sp.bar(range(len(cframe)), cframe['counts'])
            sp.xaxis.set_ticks(range(0, len(cframe)))
            sp.xaxis.set_ticklabels(cframe.index.values,
                                    rotation='vertical', ha='left')
            sp.set_xlabel("Chromosome")
            sp.set_ylabel("Frequency")
            fig.suptitle('Chromsome Counts %s' % gene)
            fig.tight_layout()
            fig.savefig('%s/chromosome_counts_%s.png' % (plot_dir, gene))

            # Plot a pie chart of the number of ERVs of each genus
            f = plt.figure()
            counts = generaL.apply(lambda x: sum(
                parsed_output['genus'] == x))
            tot = sum(counts)
            percs = counts / float(tot)
            percs = np.round(percs, 3) * 100
            percs = percs.astype(str)
            a = f.add_subplot('111')
            a.pie(counts, colors=cols)
            a.axis('equal')
            a.legend(counts.astype(str) + " " + generaL +
                     " " + percs + "%", bbox_to_anchor=(1.3, 1.1))
            f.suptitle("genera_%s" % gene)
            f.savefig('%s/genera_%s.png' % (plot_dir, gene),
                      bbox_inches='tight')

            # For each gene and genus, plot the number of ERVs assigned to
            # each group in the makeGroups step.
            j = 0
            for genus in generaL:
                bygroup = np.unique(
                    parsed_output['group'][parsed_output['genus'] == genus],
                    return_counts=True)
                if len(bygroup[0]) > 2:
                    groupf = plt.figure()
                    sp = groupf.add_subplot('111')
                    sp.bar(range(len(bygroup[1])), bygroup[1], color=cols[j],
                           edgecolor='none')
                    sp.xaxis.set_ticks(range(0, len(bygroup[1])))
                    sp.xaxis.set_ticklabels(bygroup[0], rotation='vertical',
                                            ha='left')
                    k = 0
                    for t in bygroup[1].astype(str):
                        sp.text(k, bygroup[1][k], t, va='bottom')
                        k += 1
                    groupf.suptitle("Groups %s" % genus)
                    groupf.savefig("%s/groups_%s_%s.png" % (plot_dir,
                                                            gene, genus),
                                   bbox_inches='tight')
                    j += 1
        i += 1
    out.close()
