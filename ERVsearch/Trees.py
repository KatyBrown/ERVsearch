#!/usr/bin/env python3
'''
Functions to build and process phylogenetic trees.
'''
import os
import pandas as pd
import ete3 as ete
import subprocess
import Fasta
pd.set_option('mode.chained_assignment', None)


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


def makeGroupFasta(grouptable, fastaf, path, log):
    '''
    Builds a fasta file representing each small group of ERVs using
    the output table from the assignGroups function in pipeline_ERVs.
    '''
    # make a dictionary of all reference ERVs - this is just used
    # to find the outgroup
    allfasta = Fasta.makeFastaDict("%s/ERV_db/all_ERVS.fasta" % path)
    fasta = Fasta.makeFastaDict(fastaf)
    allgroups = set(grouptable['match'])
    db_path = "%s/phylogenies" % path

    exists = dict()
    for group in allgroups:
        log.info("Making FASTA file for group %s" % group)

        f_group = dict()

        subtab = grouptable[grouptable['match'] == group]
        new = set(subtab['name'])
        for nam in new:
            f_group['%s~' % nam] = fasta[nam]

        genus, gene = group.split("_")[-2:]
        outgroup = getOutgroup(db_path, gene, genus)
        f_group['outgroup_%s' % outgroup] = allfasta[outgroup]

        group_path = "%s/group_phylogenies/%s.fasta" % (db_path, group)
        summary_path = "%s/summary_phylogenies/%s_%s.fasta" % (db_path,
                                                               gene, genus)
        if os.path.exists(group_path):
            f_group.update(Fasta.makeFastaDict(group_path))
            groupstem = "_".join(group.split("_")[:-1])
            newnam = "%s_%s" % (gene, groupstem)
            exists[newnam] = f_group
        else:
            f_group.update(Fasta.makeFastaDict(summary_path))
            exists.setdefault("%s_%s" % (gene, genus), dict())
            exists["%s_%s" % (gene, genus)].update(f_group)

    for group in exists:
        raw_out = "group_fastas.dir/%s.fasta" % (group)
        ali_out = "group_fastas.dir/%s_A.fasta" % (group)

        out = open(raw_out, "w")
        for nam, seq in exists[group].items():
            out.write(">%s\n%s\n" % (nam, seq))
        out.close()
        align(raw_out, ali_out, group, log)


def align(infile, outfile, nam, log):
    statement = ['fftns', "--quiet", infile]

    log.info("Aligning %s with FFTNS: %s %s" % (infile,
                                                " ".join(statement),
                                                outfile))
    P = subprocess.run(statement, stdin=open(infile, "r"),
                       stdout=open(outfile, "w"))
    if P.returncode != 0:
        log.error(P.stderr)
        err = RuntimeError(
            "Error aligning %s - see log file" % (nam))
        log.error(err)
        raise err


def buildTree(infile, outfile, log):
    statement = ["FastTree", "-nt", "-gtr", "-quiet", '-quote',
                 infile]
    log.info("Building phylogeny of %s: %s" % (infile, " ".join(statement)))
    P = subprocess.run(statement, stdout=open(outfile, "w"),
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr.decode())
        err = RuntimeError("Error running FastTree on %s - see log file" % (
                            infile))
        log.error(err)
        raise err


def getGroupSizes(group_names):
    gD = dict()
    for group in group_names:
        if "~" in group:
            size = len(open("group_lists.dir/%s.txt" % group[:-1]).readlines())
            gD[group] = size
    return (gD)


def drawTree(tree, outfile, maincolour, highlightcolour,
             outgroupcolour, dpi, sizenodes=False):
    '''
    Create an image file of a phylogenetic tree
    '''
    T = ete.Tree(tree, quoted_node_names=True)
    if sizenodes:
        sizeD = getGroupSizes(T.get_leaf_names())
        scale = 20 / max(sizeD.values())
    for item in T.get_leaves():
        if "outgroup" in item.name:
            T.set_outgroup(item)
            col = outgroupcolour
            text = item.name
        elif "~" in item.name:
            col = highlightcolour
            if sizenodes:
                countn = round(sizeD[item.name] * scale, 0)
                text = "%s (%i sequences)" % (item.name, sizeD[item.name])
                item.add_face(ete.CircleFace(radius=countn, color=col),
                              column=0)
            else:
                text = item.name
        else:
            col = maincolour
            text = item.name
        TF = ete.TextFace(text)
        item.add_face(TF, column=1)
        TF.fgcolor = col

    NS = ete.NodeStyle()
    NS['size'] = 0
    for node in T.traverse():
        if not node.is_leaf():
            f = ete.TextFace(node.support)
            node.add_face(f, column=0, position="branch-top")

    TS = ete.TreeStyle()
    TS.show_leaf_name = False
    T.render(outfile, tree_style=TS, dpi=int(dpi))


def monophyleticERVgroups(tree, allfasta):
    '''
    Based on a phylogenetic tree, finds monophyletic groups containing only
    novel sequences (defined by containing "~".
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
            elif "outgroup" not in leaf and leaf not in done:
                F[leaf] = allfasta[leaf]
        if count == len(leafnames):
            done = done | set(leafnames)
            sets.append(set(leafnames))
    return (sets, F)


def makeRepFastas(fastas, trees, path, outfiles, log):
    '''
    Takes a single representative for each monophyletic group of ERVs
    for monophyletic groups and builds a summary tree for each genus and gene.

    The summary tree has the summary_phylogenies sequences for the appropriate
    gene and genus, more specific reference sequences for groups where ERVs
    were identifed and a representative sequence for each new monophyletic
    ERV cluster found.
    '''
    # store sequences
    seqD = dict()
    # store number of groups in each tree (for naming)
    k = 1
    # store representaive sequences used
    repD = dict()
    phylopath = "%s/phylogenies" % path
    bn = os.path.basename(fastas[0])
    gene = bn[0:3]
    genus = os.path.splitext(bn)[0].split("_")[-1]
    allfasta = Fasta.makeFastaDict("%s/ERV_db/all_ERVS.fasta" % path)

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
