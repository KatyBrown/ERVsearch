# This file provides functions used by the pipeline.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
import ete3 as ete
import Bio.Seq as Seq
import subprocess
import copy
pd.set_option('mode.chained_assignment', None)

cols = ['crimson', 'mediumspringgreen', 'deepskyblue',
        'goldenrod', 'deeppink', 'mediumpurple', 'orangered']
genera = ['gamma', 'beta', 'spuma', 'epsilon', 'alpha', 'lenti', 'delta']
coldict = dict(zip(genera, cols))


def splitChroms(infile, log):
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
            err =  RuntimeError("""Not all chromosomes in keepchroms.txt are found in the input fasta file, the following are missing \n%s""" % ("\n".join(list(keepchroms - allchroms))))
            log.error(err)
            raise(err)
    else:
        keepchroms = allchroms

    log.info("Splitting input genome into %i chromosomes" % len(keepchroms))

    for chrom in keepchroms:
        # output a fasta file for each chromosome using samtools faidx
        statement = ["samtools", "faidx", infile, chrom]
        log.info("Processing chromosome %s: %s" % (chrom, " ".join(statement)))
        subprocess.run(statement, stdout=open(
            "host_chromosomes.dir/%s.fasta" % chrom, "w"))
        statement = ["samtools", "faidx",
                     "host_chromosomes.dir/%s.fasta" % chrom]
        log.info("Indexing chromosome %s: %s" % (chrom, " ".join(statement)))
        subprocess.run(statement)


def runExonerate(fasta, chrom, outf, log):
    '''
    Runs Exonerate protein2dna.
    Query - FASTA file containing retrovirus amino acid sequences.

    Target - chromosomes in the host_chromosomes directory.
    Settings have been optimised for time and ERV detection.
    '''

    statement = ['exonerate', '--model', 'protein2dna',
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


def getORFS(fasta, outfile_nt, outfile_orf, min_orf_len):
    '''
    Uses the EMBOSS 'getorf' function to identify open reading frames in
    the newly identified ERV regions.
    Parses the output of this software and merges it into the output dataframe
    to show the length, position and sequence of the longest ORF.
    '''

    FD = makeFastaDict(fasta)
    out_nt = open(outfile_nt, "w")
    out_orf = open(outfile_orf, "w")
    for nam, seq in FD.items():
        F = Seq.Seq(seq)
        R = F.reverse_complement()
        typs = []
        D_s = dict()
        D_t = dict()
        all_orfs = []
        for i in np.arange(0, 3):
            F_s = F[i:]
            R_s = F[i:]
            F_t = F[i:].translate()
            R_t = R[i:].translate()
            F_t_L = F_t.split("*")
            R_t_L = R_t.split("*")
            for orf in F_t_L:
                if len(orf) >= min_orf_len:
                    pos = (F_t.find(orf) * 3) + i
                    length = len(orf) * 3
                    orig_seq = F[pos: pos+length]
                    assert orig_seq.translate() == orf
                    out_orf.write(">%s_%s_%s_F+%s\n%s\n" % (nam, pos,
                                                            pos+length,
                                                            i+1, orf))
                    out_nt.write(">%s_%s_%s_F+%s\n%s\n" % (nam, pos,
                                                           pos+length,
                                                           i+1, orig_seq))
            for orf in R_t_L:
                if len(orf) >= min_orf_len:
                    pos = (R_t.find(orf) * 3) + i
                    length = len(orf) * 3
                    orig_seq = R[pos: pos+length]
                    assert orig_seq.translate() == orf
                    out_orf.write(">%s_%s_%s_F-%s\n%s\n" % (nam, pos,
                                                            pos+length,
                                                            i+1, orf))
                    out_nt.write(">%s_%s_%s_F-%s\n%s\n" % (nam, pos,
                                                           pos+length,
                                                           i+1, orig_seq))
    out_nt.close()
    out_orf.close()


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
    #  Sample groups with > 40 seqs
   # if len(groupseqs) > 40:
   #     ints = np.random.random_integers(0, len(groupseqs)-1, 20)
   #     indices = groupseqs.index.values[ints]
   #     idlist = groupseqs['id'][indices]
   # else:
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
