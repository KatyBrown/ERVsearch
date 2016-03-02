# This file provides functions used by the pipeline.

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
import re
import ete2 as ete


cols = ['crimson', 'mediumspringgreen', 'deepskyblue',
        'goldenrod', 'deeppink', 'mediumpurple', 'orangered']
genera = ['gamma', 'beta', 'spuma', 'epsilon', 'alpha', 'lenti', 'delta']
coldict = dict(zip(genera, cols))


def getParameters(ini_file):

    '''
    Reads the pipeline.ini file and generates the PARAMS dictionary of
    user-defined parameters
    '''

    D = dict()
    for line in open(ini_file).readlines():
        line = line.strip().replace(" ", "")
        if len(line) != 0:
            if line[0] != "#":
                line = line.split("=")
                D[line[0]] = line[1]
    return D


def makeChroms(infile, n_chroms):

    '''
    Takes an infile assembled into contigs or scaffolds and generates a
    set of n_chroms "chromosomes" of approximately equal length, to allow
    software designed for fully assembled genomes to run efficiently on
    partially assembled data.
    Generates chromosome_map.tsv showing the location of each scaffold and
    the host_chromosomes directory containing a fasta file representing each
    chromosome.
    '''

    mapout = open("chromosome_map.tsv", "w")

    # Count the bases in the genome
    genome_length = 0
    with open(infile) as input:
        for line in input:
            if line[0] != ">":
                genome_length += len(line.strip())

    #  Determine the approx number of bases for each chromosome
    chrom_length = genome_length / int(n_chroms)

    # Build the chromosome fasta files
    chrom_number = 1
    chrom_pos = 0
    p_chrom_pos = 0
    x = 0
    outchrom = open("host_chromosomes/chrom_%i.fa" % chrom_number, "w")
    outchrom.write(">chrom_1\n")
    contig_name = ()
    with open(infile) as input:
        for line in input:
            if line[0] == ">":
                line = line.strip()

                if x != 0:
                    # write locations to the mapout file
                    mapout.write(
                        "%s\t%i\t%i\t%i\n" % (
                            contig_name, chrom_number, p_chrom_pos, chrom_pos))
                    p_chrom_pos = chrom_pos

                # append scaffolds until the position is greater than
                # chrom_length, then start a new scaffold
                if chrom_pos >= chrom_length:
                    outchrom.close()
                    chrom_number += 1
                    outchrom = open(
                        "host_chromosomes/chrom_%i.fa" % chrom_number, "w")
                    outchrom.write(">chrom_%i\n" % chrom_number)

                    chrom_pos = 0
                    p_chrom_pos = 0
                contig_name = line.replace(">", "")
                x += 1
            else:
                outchrom.write(line)
                chrom_pos += len(line.strip())

    # write locations to the mapout file for final line
    mapout.write("%s\t%i\t%i\t%i\n" % (contig_name, chrom_number,
                                       p_chrom_pos, chrom_pos))

    outchrom.close()
    mapout.close()


def splitChroms(infile):

    '''
    For genomes assembled into chromosomes, split the input fasta file into
    one fasta file per chromosome.
    If a keep_chroms.tsv configuration file is provided, only keep the
    chromosomes listed in this file.
    '''

    outfile = ()
    if os.path.exists("keep_chroms.tsv"):
        keepchroms = [L.strip()
                      for L in open(
                      "keep_chroms.tsv").readlines()]
    else:
        keepchroms = None

    k = 0
    x = 0

    #  Generate a file each time a line starts with ">" and is in keepchroms
    #  or all lines starting with ">" if no keepchroms file is provided
    with open(infile) as input:
        for line in input:
            if line[0] == ">":
                chromname = line.replace(">", "").strip()
                if x != 0:
                    outfile.close()
                    k = 0
                if (keepchroms is None or chromname in keepchroms):
                        k = 1
                        s = chromname.replace("chr", "")
                        s = s.replace("chrom", "")
                        s = s.replace("_", "-")
                        cname = "chrom_%s" % s
                        outfile = open(
                            "host_chromosomes/%s.fa" % cname, "w")
                        outfile.write(">%s\n" % cname)
                        x += 1
            elif k == 1:
                outfile.write(line)
    outfile.close()


def runExonerate(fasta, chrom, out, path):

    '''
    Runs Exonerate protein2dna.
    Query - retrovirus sequences in the "sequencedir" directory provided in
    pipeline.ini.  Must be fasta files of ERV sequences ending with .fa.
    Recommended to use the provided ERV sequence files.
    Target - chromosomes in the host_chromosomes directory.
    Settings have been optimised for time and ERV detection.
    '''

    statement = """
    %s \
    --model protein2dna \
    --showalignment F \
    --seedrepeat 1 \
    --showvulgar T \
    --query %s \
    --target %s \
    > %s""" % (path, fasta, chrom, out)
    os.system(statement)


def filterExonerate(out, min_hit_length):

    '''
    Filters hits shorter than the min_hit_length specified in pipeline.ini
    from the Exonerate output.  Parses the columns in the Exonerate output
    table and generates a pandas dataframe of this information - the
    results output table.
    '''

    # Parsed column names.  Details are from the final column of the
    # Exonerate showvulgar output, as listed here
    # http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual

    cnames = ["query_id", "query_start", "query_end", "query_strand",
              "target_id", "target_start", "target_end",
              "target_strand", "score", "details", "length"]
    smalldf = pd.DataFrame(columns=cnames)
    j = 0

    # Read the input file
    with open(out) as infile:
        for line in infile:
            line = line.strip().split(": ")
            if line[0] == "vulgar":
                segs = line[1].split(" ")
                for i in range(9):
                    smalldf.ix[str(j), cnames[i]] = segs[i]
                rest = "|".join(segs[9:])
                smalldf.ix[str(j), cnames[9]] = rest
                j += 1

    # Find the hit lengths (allowing for strand)
    smalldf['target_start'] = smalldf['target_start'].astype(int)
    smalldf['target_end'] = smalldf['target_end'].astype(int)

    smalldf['length'] = abs(smalldf['target_end'] -
                            smalldf['target_start'])

    plus = smalldf[smalldf['target_start'] <= smalldf['target_end']]
    minus = smalldf[smalldf['target_start'] > smalldf['target_end']]

    tempstart = minus['target_start']
    tempend = minus['target_end']

    minus = minus.drop(['target_start', 'target_end'], 1)
    minus['target_start'] = tempend
    minus['target_end'] = tempstart

    smalldf = plus.append(minus)
    smalldf = smalldf.sort(['target_start'])

    # Filter hits which are too small or contain introns
    smalldf = smalldf[((smalldf['length'] > min_hit_length) &
                       (smalldf.details.str.contains("I") == False))]

    smalldf['target_start'] = smalldf['target_start'].astype('object')
    smalldf['target_end'] = smalldf['target_end'].astype('object')

    return smalldf[cnames]


def makeFastaDict(multifasta):

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


def getORFS(matches, fasta, outfiles, emboss):

    '''
    Uses the EMBOSS 'getorf' function to identify open reading frames in
    the newly identified ERV regions.
    Parses the output of this software and merges it into the output dataframe
    to show the length, position and sequence of the longest ORF.
    '''

    data = pd.read_csv(matches, sep="\t", index_col=0)

    #  Run EMBOSS getorf
    if len(data) != 0:
        fasout = outfiles[1]
        os.system("""sed 's/\:/_/g' %s |
                     %s/getorf -filter -outseq %s""" % (fasta, emboss, fasout))
        fasta = makeFastaDict(fasout)

        # Find and record the starts, ends, lengths and seqeunces of the ORFs
        # in the getorf output.
        res = []
        uniqs = []
        for key in fasta.keys():
            segs = re.split("_|\[|\]|\s+", key.strip())
            nam = "chrom_%s:%s" % (segs[1], segs[2])
            uniq = "_".join(segs[0:5])
            uniqs.append(uniq)
            div = segs[2].split("-")
            ss, se = int(div[0]), int(div[1])
            length = abs(ss - se)
            all = [nam, ss, se, length, uniq, fasta[key]]
            res.append(all)

        # Make a temporary dataframe of this information and sort it by
        # length for each ERV then take the longest only.
        R = pd.DataFrame(res, columns=[
            'id', 'start', 'end', 'length', 'uniq', 'seq'], index=uniqs)
        R = R.sort_values(['id', 'length'])
        R = R.groupby('id').last()
        R = R.drop('uniq', 1)
        R.columns = ['orf_start', 'orf_end', 'orf_len', 'orf_seq']

        # Merge this into the results dataframe.
        data = data.merge(R, left_index=True, right_index=True)
        data.to_csv(outfiles[0], sep="\t")
    else:
        os.system("touch %s; touch %s" % (outfiles[0], outfiles[1]))


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
    if len(groupseqs) > 40:
        ints = np.random.random_integers(0, len(groupseqs)-1, 20)
        indices = groupseqs.index.values[ints]
        idlist = groupseqs['id'][indices]
    else:
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
        os.system("%s/fftns --quiet %s > %s"
                  % (mafftpath, outf, outf.replace(".fasta", "_ali.fasta")))

        # Build Trees
        os.system("%s/FastTree -nt -gtr -quiet %s > %s" % (fasttreepath,
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
    novel sequences (defined by the "chrom" prefix.
    '''
    p = ete.Tree(tree)
    sets = []
    all = []

    # Runs through all the subtrees of the specified tree and determines
    # if they contain any non-novel groups.
    for x in p.traverse():
        leafnames = x.get_leaf_names()
        count = sum([1 if leaf.startswith("chrom_") else 0
                     for leaf in leafnames])
        all.append([leaf if leaf.startswith("chrom_") else ''
                    for leaf in leafnames])
        if count == len(leafnames):
            sets.append(set(leafnames))
    all = set([b for a in all for b in a])
    all.remove('')

    # sort the groups by desending size (by number of elements)
    pairs = [(s, len(s)) for s in sets]
    sortedgroups = [a[0]
                    for a in sorted(pairs, key=lambda tup: tup[1])[::-1]]
    i = 0
    # keep only groups which are not a subset of another group
    while i < len(sortedgroups):
        sortedgroups = remove_subsets(sortedgroups[i], sortedgroups)
        i += 1

    return sortedgroups


def correctGroupSize(groupsizes, tree, genetab):
    '''
    Corrects the group size for sequences representing a group if the original
    tree was based on a sample of sequences, rather than all the sequences
    in a group
    '''
    total = sum(groupsizes)
    treenam = tree.replace(".tre", "").replace("trees/", "")

    allmatches1 = set(genetab['id'][genetab['group'] == treenam])
    allmatches2 = set(genetab['id'][genetab['match'] == treenam])
    allmatches = allmatches1 | allmatches2

    if len(allmatches) == total:
        tt = "full"
        return groupsizes, tt

    elif len(allmatches) >= 40:
        tt = "approx"
        n_represented = int(round((float(len(allmatches)) / 20), 0))
        cg = np.array(groupsizes) * n_represented
        return cg, tt


def makeRepFastas(infiles, genetab, wd, pp,
                  mafftpath, fasttreepath, seqdir, outfile):
    '''
    Using the Fasta files from makePhyloFastas, takes a single representative
    for monophyletic groups and builds a summary tree for each genus and gene.
    Generates a phylogeny and an image of each phylogeny with node labels sized
    according to the size of the group.
    '''
    allfasta = makeFastaDict("%s/all_ERVS.fasta" % seqdir)
    repcounts = dict()
    repseqs = dict()
    repgroups = dict()

    # Choose a representative for each monophyletic cluster of novel ERV
    # regions in the input trees
    for tree in infiles:
        fasta = tree.replace("tree", "fasta")
        fasta = fasta.replace(".tre", ".fasta")
        fasta = makeFastaDict(fasta)

        groups = monophyleticERVgroups(tree)
        groupsizes = [len(group) for group in groups]

        # correct group sizes - for very large groups only a subset is in
        # the input tree
        corrected_groupsizes, tt = correctGroupSize(groupsizes, tree, genetab)

        # Record the representative member of each group, its sequence and
        # the corrected group size.  If the input trees contain all
        # sequences, rather than subsets of sequences, make an outfile
        # with an ID number containing all the sequence names in the group

        j = 0
        for group in groups:
            L = list(group)
            rep = L[0]
            groupsizes.append(len(L))
            rep = L[0]
            repcounts[rep] = corrected_groupsizes[j]
            repseqs[rep] = fasta[rep]

            if tt == "full":
                groupid = str(id(group))[-6:]
                gout = open("group_lists/%s.txt" % groupid, "w")
                for nam in group:
                    if nam == rep:
                        gout.write("%s**\n" % nam)
                    else:
                        gout.write("%s\n" % nam)
                gout.close()
                repgroups[rep] = groupid
            else:
                repgroups[rep] = ""
            j += 1

    # Build a fasta file with the group representatives and the representative
    # reference trees
    genus, gene = outfile.split("_")[-2:]
    gene = gene.replace(".tre", "")
    genus = genus.replace("trees/", "")
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
    os.system("%s/fftns --quiet %s > %s"
              % (mafftpath, fastaout,
                 fastaout.replace(".fasta", "_ali.fasta")))

    # Build a tree of the fasta file
    os.system("%s/FastTree -nt -gtr -quiet %s > %s"
              % (fasttreepath, aliout, outfile))

    # Build an image of the fasta file
    treepng = outfile.replace(".tre", ".png")
    T = ete.Tree(outfile)
    m = max(repcounts.values())
    scale = 40 / float(m)

    for item in T.traverse():
        if item.name[0:3] == "chr":
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
            out.write("Mean Length %s = %f\n" % (gene,
                                                 np.mean(
                                                     parsed_output['length'])))
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
