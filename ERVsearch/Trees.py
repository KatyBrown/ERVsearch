#!/usr/bin/env python3
'''
Functions to build and process phylogenetic trees.
'''
import os
import pandas as pd
import ete3 as ete
import subprocess
try:
    import ERVsearch.Fasta as Fasta
    import ERVsearch.Errors as Errors
except ImportError:
    import Fasta
    import Errors
pd.set_option('mode.chained_assignment', None)


def getOutgroup(phylopath, gene, genus):
    '''
    Uses the provided reference file to identify an appropriate outgroup
    for a particular gene and genus of ERVs.
    Returns the outgroup sequence name.
    '''
    # The outgroups file is just a list of sequences to use as outgroups
    # for each gene and genus
    outgdict = [tuple(line.strip().split("\t"))
                for line in open(
                        "%s/outgroups.tsv" % phylopath).readlines()]
    outgdict = dict(outgdict)
    genegenus = "%s_%s" % (genus, gene)
    outg = outgdict[genegenus]
    return outg


def makeGroupFasta(grouptable, fastaf, path, log):
    '''
    Builds a fasta file representing each small group of ERVs using
    the output table from the assignGroups function in the main ERVsearch
    pipeline, an existing FASTA file of sequences in the group of interest
    and the newly identified ERV-like sequences.
    '''
    # make a dictionary of all reference ERVs - this is just used
    # to find the outgroup sequence
    allfasta = Fasta.makeFastaDict("%s/ERV_db/all_ERVs_nt.fasta" % path)

    # make a dictionary of the new ERVs
    fasta = Fasta.makeFastaDict(fastaf)

    # find all the groups any sequences have been assigned to
    allgroups = set(grouptable['group'])
    db_path = "%s/phylogenies" % path

    exists = dict()
    # for each group any sequence has been assigned to
    for group in allgroups:
        log.info("Making FASTA file for group %s" % group)

        f_group = dict()

        # get the names of all the new sequences in the group
        subtab = grouptable[grouptable['group'] == group]
        new = set(subtab['name'])

        # put the new sequences into the dictionary
        for nam in new:
            # Add "~" to the new sequence names so they are easy to find
            f_group['%s~' % nam] = fasta[nam]

        # get the gene and genus
        genus, gene = group.split("_")[-2:]

        # find an appropriate outgroup
        outgroup = getOutgroup(db_path, gene, genus)

        # put the outgroup in the dictionary
        f_group['outgroup_%s' % outgroup] = allfasta[outgroup]

        # paths to the reference fasta files
        group_path = "%s/group_phylogenies/%s.fasta" % (db_path, group)
        summary_path = "%s/summary_phylogenies/%s_%s.fasta" % (db_path,
                                                               gene, genus)

        # if there is a reference file for this group, put the sequences in
        # dictionary
        if os.path.exists(group_path):
            f_group.update(Fasta.makeFastaDict(group_path))
            groupstem = "_".join(group.split("_")[:-1])
            newnam = "%s_%s" % (gene, groupstem)
            # store the dictionary for the group
            exists[newnam] = f_group
        # otherwise use the reference file for the genus
        else:
            f_group.update(Fasta.makeFastaDict(summary_path))
            exists.setdefault("%s_%s" % (gene, genus), dict())
            # store the dictionary for the group
            exists["%s_%s" % (gene, genus)].update(f_group)

    for group in exists:
        # unaligned output file
        raw_out = "group_fastas.dir/%s.fasta" % (group)
        # aligned output file
        ali_out = "group_fastas.dir/%s_A.fasta" % (group)

        log.info("Writing output FASTA for group %s" % group)

        # write the sequences to the output file
        out = open(raw_out, "w")
        for nam, seq in exists[group].items():
            out.write(">%s\n%s\n" % (nam, seq))
        out.close()
        # align them
        align(raw_out, ali_out, group, log)


def align(infile, outfile, nam, log):
    '''
    Aligns a FASTA file using the default settings of the MAFFT fftns
    algorithm (https://mafft.cbrc.jp/alignment/software/manual/manual.html)
    '''
    # Build the MAFFT statement
    statement = ['fftns', "--quiet", infile]
    log.info("Aligning %s with FFTNS: %s %s" % (infile,
                                                " ".join(statement),
                                                outfile))

    # Run MAFFT
    P = subprocess.run(statement, stdin=open(infile, "r"),
                       stdout=open(outfile, "w"))
    if P.returncode != 0:
        log.error(P.stderr)
        Errors.raiseError(Errors.SoftwareError, infile, "mafft", log=log)


def buildTree(infile, outfile, log):
    '''
    Builds a phylogenetic tree using the FastTree2 algorithm
    (http://www.microbesonline.org/fasttree) for a nucleotide
    input file, using the default settings plus the GTR model.
    '''
    # Build the FastTree statement
    # Quote is used because of the punctuation in the names
    statement = ["FastTree", "-nt", "-gtr", "-quiet", '-quote',
                 infile]
    log.info("Building phylogeny of %s: %s" % (infile, " ".join(statement)))

    # Run FastTree
    P = subprocess.run(statement, stdout=open(outfile, "w"),
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr.decode())
        Errors.raiseError(Errors.SoftwareError, infile, "FastTree", log=log)


def getGroupSizes(group_names):
    '''
    Determines the sizes of the groups used to scale the nodes in the
    phylogenetic tree image files.

    The sizes are found by counting how many sequences are in the relevant
    file in group_lists.dir.

    A dictionary is returned where keys are group names and values are
    relative sizes.
    '''
    gD = dict()
    for group in group_names:
        # used to mark newly identified sequences
        if "~" in group:
            size = len(open("group_lists.dir/%s.txt" % group[:-1]).readlines())
            gD[group] = size
    return (gD)


def drawTree(tree, outfile, maincolour, highlightcolour,
             outgroupcolour, dpi, sizenodes=False):
    '''
    Create an image file of a phylogenetic tree.
    Leaf names are in maincolour unless the leaf name contains a "~" in which
    case they are in highlightcolour. The outgroup (which has outgroup in the
    name) is labelled in outgroupcolour.
    The output file dpi is determined using the dpi parameter.

    If sizenodes is True, the nodes are scaled by the numer of sequences they
    represent. If it is False all the nodes are the same size.
    '''
    # Read the tree into ete3
    T = ete.Tree(tree, quoted_node_names=True)
    if sizenodes:
        # The node sizes are established by counting how many sequences
        # are in the group in the group_lists.dir directory.
        sizeD = getGroupSizes(T.get_leaf_names())
        # This scaling seems to make the nodes an appropriate size to fit
        # on the page
        scale = 20 / max(sizeD.values())

    for item in T.get_leaves():
        # Look for the outgroup and set it as the outgroup
        if "outgroup" in item.name:
            T.set_outgroup(item)
            # set the outgroup colour
            col = outgroupcolour
            text = item.name
        elif "~" in item.name:
            col = highlightcolour
            if sizenodes:
                # scale the nodes if needed
                countn = round(sizeD[item.name] * scale, 0)
                # Write the sequence name and the number of sequences it
                # is representing
                text = "%s (%i sequences)" % (item.name, sizeD[item.name])
                item.add_face(ete.CircleFace(radius=countn, color=col),
                              column=0)
            else:
                # Otherwise just put the name
                text = item.name
        else:
            col = maincolour
            text = item.name

        # Add the text to the node
        TF = ete.TextFace(text)
        item.add_face(TF, column=1)
        TF.fgcolor = col

    # This stops the default behaviour in ete3 to have every node (including
    # internal nodes) to be marked with a dot.
    NS = ete.NodeStyle()
    NS['size'] = 0
    for node in T.traverse():
        if not node.is_leaf():
            f = ete.TextFace(node.support)
            # This sets the position of the branch support
            node.add_face(f, column=0, position="branch-top")
        node.set_style(NS)

    TS = ete.TreeStyle()
    # Hide the original leaf names because we have reassigned them
    TS.show_leaf_name = False

    # Make the output file
    T.render(outfile, tree_style=TS, dpi=int(dpi))


def monophyleticERVgroups(tree, allfasta):
    '''
    Based on a phylogenetic tree, finds monophyletic groups containing only
    novel sequences (defined by containing "~").

    Only the largest monophyletic group containing each ERV is returned.
    '''
    sets = []
    done = set()
    F = dict()
    # Runs through all the subtrees of the specified tree and determines
    # if they contain any non-novel groups.
    for x in tree.traverse():
        leafnames = x.get_leaf_names()
        count = 0
        for leaf in leafnames:
            if "~" in leaf and leaf not in done:
                count += 1
            # also keep the names of the reference sequences
            elif "outgroup" not in leaf and leaf not in done:
                F[leaf] = allfasta[leaf]
        if count == len(leafnames):
            done = done | set(leafnames)
            sets.append(set(leafnames))
    return (sets, F)


def makeRepFastas(fastas, trees, path, outfiles, log):
    '''
    Takes matched pairs of FASTA files and trees and choses a single
    representative sequence (the longest) for each monophyletic group
    of newly identified ERVs.

    A summary fasta file is then generated containing the summary_phylogenies
    sequences for the appropriate gene and genus, more specific reference
    sequences for groups where ERVs were identifed and a representative
    sequence for each new monophyletic ERV cluster found.
    '''
    # store sequences
    seqD = dict()
    # store number of groups in each tree (for naming)
    k = 1
    # store representaive sequences used
    repD = dict()

    # path to the phylogenies directory
    phylopath = "%s/phylogenies" % path
    bn = os.path.basename(fastas[0])
    gene = bn[0:3]
    genus = os.path.splitext(bn)[0].split("_")[-1]

    # put all ERV sequences into a dictionary - this is just used
    # to get the outgroup sequence
    allfasta = Fasta.makeFastaDict("%s/ERV_db/all_ERVs_nt.fasta" % path)

    for fasta, tree in zip(fastas, trees):
        log.info("Incorporating %s into summary phylogeny" % tree)

        # convert the fasta to a dictionary and read the tree with ete3
        F = Fasta.makeFastaDict(fasta)
        T = ete.Tree(tree, quoted_node_names=True)
        # get the local reference sequences and the monophyletic ERV groups
        groups, refs = monophyleticERVgroups(T, allfasta)

        log.info("Identified %s monophyletic ERV groups in %s" % (len(groups),
                                                                  tree))

        # store the local reference sequences for this group
        seqD.update(refs)

        for group in groups:
            seqs = [F[nam] for nam in group]
            # sort from longest to smallest
            Z = sorted(zip(group, seqs), key=lambda x: len(x[1]),
                       reverse=True)
            # use the longest as the representative
            rep = list(Z)[0][0]
            ID = "%s_%s_%i~" % (gene, genus, k)

            # write the sequence names and sequences to file
            gout = open("group_lists.dir/%s.txt" % ID[:-1], "w")
            gfasta = open("group_lists.dir/%s.fasta" % ID[:-1], "w")
            for nam in group:
                gfasta.write(">%s\n%s\n" % (nam, F[nam]))
                if nam == rep:
                    gout.write("%s**\n" % nam)
                else:
                    gout.write("%s\n" % nam)
            gout.close()
            gfasta.close()

            # update the results
            repD[ID] = rep
            seqD[ID] = F[rep]
            k += 1

    # write the representative sequences to file
    reps = open("group_lists.dir/%s_%s_representative_sequences.txt" % (
        gene, genus), "w")
    for nam, rep in repD.items():
        reps.write("%s\t%s\n" % (nam, rep))
    reps.close()
    # Build a fasta file with the group representatives and the representative
    # reference trees

    # read in the summary reference sequences
    allD = Fasta.makeFastaDict("%s/summary_phylogenies/%s_%s.fasta" % (
        phylopath, gene, genus))
    allD.update(seqD)
    # find the outgroup
    outg = getOutgroup(phylopath, gene, genus)
    allD["outgroup_%s" % outg] = allfasta[outg]

    # write the output to file
    out = open(outfiles[0], "w")
    for nam, seq in allD.items():
        out.write(">%s\n%s\n" % (nam, seq))
    out.close()

    # align
    align(outfiles[0], outfiles[1], "%s_%s" % (gene, genus), log)
