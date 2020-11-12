# This is the main ruffus pipeline for ERVsearch.

from ruffus import follows, split, transform, mkdir, formatter, originate
from ruffus import suffix, collate, regex, merge, add_inputs
from ruffus.combinatorics import product

import PipelineERVs as PipelineERVs
import os
import pandas as pd
import numpy as np
import ruffus.cmdline as cmdline
import logging
import subprocess
import configparser
import pybedtools

global PARAMS
global log


parser = cmdline.get_argparse(description='Pipeline ERVs')
options = parser.parse_args()
PARAMS = configparser.ConfigParser()
PARAMS.read("pipeline.ini")

outstem = PARAMS['output']['outfile_stem']
# Set up logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
logfile = "%s_log.txt" % outstem
handler = logging.FileHandler(logfile)
handler.setLevel(logging.INFO)
formats = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formats)
log.addHandler(handler)

genes = []
PARAMS['gene'] = {}

usecustom = PARAMS['database']['use_custom_db'] == "True"

for gene in ['gag', 'pol', 'env']:
    if PARAMS['database'][gene] != "None" or not usecustom:
        if usecustom:
            genedb = PARAMS['database'][gene]
        else:
            genedb = "%s/ERV_db/%ss.fa" % (
                PARAMS['database']['path_to_ERVsearch'], gene)
        if not os.path.exists(genedb):
            err = RuntimeError(
                "Database fasta file %s does not exist" % genedb)
            log.error(err)
            raise err
        PARAMS['gene'][gene] = genedb
        genes.append(gene)

PARAMS.write(open("%s_parameters.txt" % outstem, "w"))


@follows(mkdir("host_chromosomes.dir"))
@split(PARAMS['input']['genome'], r"host_chromosomes.dir/*.fasta")
def genomeToChroms(infile, outfiles):
    '''
    Splits the host genome provided by the user into one fasta file for each
    chromosome, scaffold or contig.

    A temporary unzipped copy of zipped and gzipped fasta files will be
    created.

    This function generates a series of fasta files which are stored in the
    host_chromosomes.dir directory.
    '''
    assert os.path.exists(infile), "Input file %s not found" % infile
    if infile.endswith(".gz"):
        # unzip gzipped files
        log.info("Unzipping gzipped input file %s" % infile)
        stem = os.path.basename(infile).split(".gz")[0]
        statement = ["gunzip", "-c", infile]
        new_infile = stem
        log.info("Running statement: %s" % " ".join(statement))
        subprocess.run(statement, stdout=open(new_infile, "w"))
    elif infile.endswith(".zip"):
        # unzip zipped files
        log.info("Unzipping zipped input file %s" % infile)
        stem = os.path.basename(infile).split(".zip")[0]
        statement = ["unzip", "-p", infile]
        new_infile = stem
        log.info("Running statement: %s" % " ".join(statement))
        subprocess.run(statement, stdout=open(new_infile, "w"))
    else:
        log.info("Linking to input file %s" % infile)
        stem = os.path.basename(infile)
        statement = ["ln",  "-sf", infile]
        subprocess.run(statement)
        new_infile = stem

    statement = ["samtools", "faidx", new_infile]
    log.info("Indexing fasta file: %s" % " ".join(statement))
    s = subprocess.run(statement)
    if s.returncode != 0:
        err = RuntimeError("""Indexing the input Fasta file was not possible, please check your input file for errors""")
        log.error(err)
        raise(err)
    PipelineERVs.splitChroms(new_infile, log)
    log.info("Removing temporary input file %s" % new_infile)
    os.unlink(new_infile)
    os.unlink("%s.fai" % new_infile)


@follows(mkdir("raw_exonerate_output.dir"))
@follows(mkdir("clean_exonerate_output.dir"))
@follows(genomeToChroms)
@product(genomeToChroms, formatter(),
         [PARAMS['gene'][gene] for gene in genes],
         formatter(),
         r'raw_exonerate_output.dir/{basename[1][0]}_{basename[0][0]}.tsv')
def runExonerate(infiles, outfiles):
    '''
    Runs the protein2dna algorithm in the Exonerate software package with
    the host chromosomes in host_chromosomes as target sequences and the
    ERV_Amino_Acid_DB fasta files as query sequences.

    Output is parsed to remove sequences shorter than the "overlap" parameter
    in pipeline.ini and remove all sequences containing introns and
    combined into a tab-delimited file

    The raw output of Exonerate is stored in the raw_exonerate_output directory
    This step is carried out with low stringency as results are later filtered
    using UBLAST and Exonerate.

    The filtered output is stored in the clean_exonerate_output directory as
    gag_XXX_filtered.tsv, pol_XXX_filtered.tsv and env_XXX_filtered.tsv,
    where XXX is the chromsome or contig name.

    A bed file (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    is also generated corresponding to the regions in each filtered output
    file - these are stored in clean_exonerate_output as
    gag_XXX_filtered.bed, pol_XXX_filtered.bed, env_XXX_filtered.bed.
    '''

    log.info("Running Exonerate on %s vs %s" % (infiles[0], infiles[1]))

    PipelineERVs.runExonerate(infiles[1], infiles[0],
                              outfiles[0], log)


@transform(runExonerate, regex("raw_exonerate_output.dir/(.*).tsv"),
           [r'clean_exonerate_output.dir/\1_unfiltered.tsv',
            r'clean_exonerate_output.dir/\1_filtered.tsv',
            r'clean_exonerate_output.dir/\1.bed'])
def cleanExonerate(infile, outfiles):
    # Minimum size number of aas in a hit t
    min_hit_length = int(PARAMS['exonerate']['min_hit_length'])
    PipelineERVs.filterExonerate(infile, outfiles, min_hit_length, log)


@follows(mkdir("gene_bed_files.dir"))
@collate(cleanExonerate,
         regex("clean_exonerate_output.dir/([a-z]+)s_(.*)_unfiltered.tsv"),
         [r"gene_bed_files.dir/\1_all.bed",
          r"gene_bed_files.dir/\1_merged.bed"])
def mergeOverlaps(infiles, outfiles):
    '''
    Overlapping regions of the genome detected by Exonerate with
    similarity to the same retroviral gene are merged into
    single regions.  This is performed using bedtools merge on
    the bed files output by cleanExonerate.
    Merged bed files are stored as
    gags_merged.bed, pols_merged.bed and envs_merged.bed.
    '''
    log.info("Generating combined bed file %s" % outfiles[0])
    beds = [inf[2] for inf in infiles]
    combined = PipelineERVs.combineBeds(beds)
    log.info("%i records identified in combined bed file %s" % (
        len(combined), outfiles[0]))
    combined.to_csv(outfiles[0], sep="\t", index=None, header=None)

    merged = PipelineERVs.mergeBed(outfiles[0], PARAMS['exonerate']['overlap'],
                                   log)
    log.info("Writing merged bed file %s with %i lines" % (
        outfiles[1], len(merged)))
    merged.to_csv(outfiles[1], sep="\t", index=None, header=None)


@split(mergeOverlaps, ("gag_merged.fasta",
                       "pol_merged.fasta",
                       "env_merged.fasta"))
def makeFastas(infiles, outfile):
    '''
    Fasta files are generated containing the sequences
    of the merged regions of the genome identified using
    mergeOverlaps.
    These are extracted from the host chromosomes using Samtools.
    The output files are stored in clean_exonerate_output as
    gags.fa, pols.fa and envs.fa.
    '''
    for infile in infiles:
        if "merged" in infile:
            outfile = infile.replace("_merged.bed", "_merged.fasta")
            statement = 'bedtools getfasta -s -fi %s/%s.fa -bed %s > %s' % (
                    PARAMS['genome_directory'], PARAMS['genome'],
                    infile, outfile)
            log.info(statement)
            os.system(statement)

@transform(makeFastas, suffix("_merged.fasta"),
           "_unfiltered_exonerate_merged.fasta")
def renameFastas(infile, outfile):
    out = open(outfile, "w")
    F = PipelineERVs.makeFastaDict(infile)
    gene = infile.split("_")[0]
    for i, nam in enumerate(F):
        out.write(">%s%s_%s\n%s\n" % (gene, i+1, nam, F[nam]))
    out.close()


@transform("%s/*fa" % PARAMS["database"],
           regex("%s/(.*).fa" % PARAMS["database"]),
           r'UBLAST_db/\1.udb')
def makeUBLASTDb(infile, outfile):
    '''
    USEARCH requires an indexed database of query sequences to run.
    This function generates this database for the three
    ERV_Amino_Acid_DB fasta files.
    '''
    os.system("%s -makeudb_ublast \
               %s -output %s  -quiet" % (
                   PARAMS['path_to_usearch'], infile, outfile))


@transform(renameFastas, suffix("_unfiltered_exonerate_merged.fasta"),
           add_inputs(makeUBLASTDb),
           ("_filtered_UBLAST.tsv", "_raw_UBLAST_alignments.txt", "_filtered_UBLAST.fa"))
def runUBLASTCheck(infiles, outfiles):
    '''
    ERV regions in the fasta files generated by makeFasta
    are compared to the ERV_Amino_Acid_DB files for a second
    time, this time using USEARCH.
    This allows sequences with low similarity to known ERVs
    to be filtered out.  Similarity thresholds can be set in
    the pipeline.ini file (usearch_id, min_hit_length and usearch_coverage).
    The output file is a fasta file of sequences with high similarity
    to known retroviruses in the ERV_Amino_Acid_DB.
    These are saved in parsed_exonerate_output as
    gags_filtered.fa, pols_filtered.fa and envs_filtered.fa.
    Raw output is also saved in parsed_exonerate_output as
    gags_alignments.txt, pols_alignments.txt and envs_alignments.txt.
    '''
    fasta = infiles[0]
    fastanam = fasta.split("/")[-1].split("_")[0] + "s"
    for inf in infiles[1:]:
        if inf.split("/")[-1].replace(".udb", "") == fastanam:
            dbnam = inf
    id_thresh = PARAMS['usearch_id']
    min_coverage = PARAMS['usearch_coverage']
    os.system("""%s -ublast %s \
                 -db %s \
                 -evalue 1 -query_cov %s -id %s\
                 -top_hit_only \
                 -blast6out %s \
                 -quiet \
                 -alnout %s""" % (PARAMS['path_to_usearch'], fasta, dbnam,
                                  min_coverage, id_thresh, outfiles[0],
                                  outfiles[1]))
    L = set([line.strip().split("\t")[0]
             for line in open(outfiles[0]).readlines()])
    PipelineERVs.simpleFasta(L, infiles[0], outfiles[2])


@transform(runUBLASTCheck, suffix("_filtered_UBLAST.tsv"),
           add_inputs(PARAMS['database']),
           ["_all_matches_exonerate.tsv",
            "_best_matches_exonerate.tsv",
            "_refiltered_exonerate.fasta"])
def findBest(infiles, outfiles):
    '''
    Runs the exonerate ungapped algorithm with each ERV region
    in the fasta files generated by makeFasta as queries and the
    All_ERVs_Fasta fasta file as a target, to detect which known
    retrovirus is most similar to each newly identified ERV region.
    All_ERVs_Fasta contains nucleic acid sequences for many known
    endogenous and exogenous retroviruses
    The raw output is saved in parsed_exonerate_output as
    gags_table_matches.tsv, pols_table_matches.tsv and
    envs_table_matches.tsv.
    Regions with no significant similarity to a known retrovirus
    are filtered out.
    The most similar known retrovirus to each of the ERV regions is
    identified.
    This result is saved as gags_table_bestmatches.tsv,
    pols_table_bestmatches.tsv and envs_table_bestmatches.tsv
    in the parsed_exonerate_output directory.
    '''
    gene = infiles[0][0].split("/")[-1].split("_")[0]
    fasta = infiles[0][2]
    refs = infiles[1]
    os.system("""
    %s \
    --query %s \
    --target %s \
    --showalignment F \
    --showvulgar F \
    --ryo "%%qi\t%%ti\t%%s\n" \
    --verbose 0 \
    --score %s\
    | sort -n \
    > %s """ % (PARAMS['path_to_exonerate'],
                fasta, refs, PARAMS['exonerate_minscore'],
                outfiles[0]))
    res = pd.read_csv(outfiles[0],
                      sep="\t", header=None, names=['id', 'match', 'score'])
    res['gene'] = res['match'].str.split("_").str.get(-1)
    res = res[res['gene'] == gene]
    res = res.drop('gene', 1)
    rtab = pd.DataFrame(columns=['id', 'match', 'score'])
    ids = np.unique(res['id'].values)
    for id in ids:
        subdf = res[res['id'] == id]
        maxs = subdf[subdf['score'] == np.max(subdf['score'])].values[0]
        rtab.loc[maxs[0]] = maxs
    rtab['gene'] = rtab['match'].str.split("_").str.get(-1)
    rtab['genus'] = rtab['match'].str.split("_").str.get(-2)
    rtab.to_csv(outfiles[1], sep="\t", index=None)

    L = list(set(rtab['id']))
    PipelineERVs.simpleFasta(L, fasta, outfiles[2])


@transform(findBest, suffix("_all_matches_exonerate.tsv"),
           ["_orfs_nt.fasta", "_orfs_aa.fasta"])
def ORFs(infiles, outfiles):
    '''
    Finds the longest open reading frame in each of the ERV regions
    in the output table
    This analysis is performed using EMBOSS getorfs
    Raw getorfs output is saved in parsed_exonerate_output as
    gags_table_orfs.fa, pols_table_orfs.fa and envs_table_orfs.fa.
    The start, end, length and sequence of each ORF are added to the
    output tables and saved in parsed_exonerate_output as
    gags_table_orfs.tsv, pols_table_orfs.tsv and envs_table_orfs.tsv.
    '''
    fasta = infiles[2]
    PipelineERVs.getORFS(fasta, outfiles[0], outfiles[1],
                         int(PARAMS['min_orf_len']))


@transform(ORFs, suffix("_orfs_nt.fasta"),
           add_inputs(makeUBLASTDb),
           ("_filtered_UBLAST_ORFs.tsv",
            "_raw_UBLAST_alignments_ORFs.txt",
            "_filtered_UBLAST_ORFs_aa.fasta",
            "_filtered_UBLAST_ORFs_nt.fasta"))
def checkORFsUBLAST(infiles, outfiles):
    fasta = infiles[0][1]
    fastanam = fasta.split("/")[-1].split("_")[0] + "s"
    for inf in infiles[1:]:
        if inf.split("/")[-1].replace(".udb", "") == fastanam:
            dbnam = inf
    id_thresh = PARAMS['usearch_id']
    min_coverage = PARAMS['usearch_coverage']
    os.system("""%s -ublast %s \
                 -db %s \
                 -evalue 1 -query_cov %s -id %s\
                 -top_hit_only \
                 -blast6out %s \
                 -quiet \
                 -alnout %s""" % (PARAMS['path_to_usearch'], fasta, dbnam,
                                  min_coverage, id_thresh, outfiles[0],
                                  outfiles[1]))
    L = set([line.strip().split("\t")[0]
             for line in open(outfiles[0]).readlines()])
    PipelineERVs.simpleFasta(L, infiles[0][0], outfiles[3])
    PipelineERVs.simpleFasta(L, infiles[0][1], outfiles[2])


@transform(checkORFsUBLAST, suffix("_filtered_UBLAST_ORFs.tsv"),
           add_inputs(findBest),
           "_groups.tsv")
def makeGroups(infiles, outfile):
    '''
    The retroviruses in All_ERVs_Fasta have been classified
    into groups based on sequence similarity.
    Each group is named after a single representative ERV.
    The newly identified ERV regions are classified into the
    same groups based on the output of the findBest function.
    Each region is assigned to the same group as the retrovirus
    sequence it was found to be most similar to by findBest.
    The assigned group is added to the output tables, these are
    saved as gags_table_groups.tsv, pols_table_groups.tsv and
    envs_table_groups.tsv.
    '''
    ORFs = infiles[0][2]
    ORF_fasta = PipelineERVs.makeFastaDict(ORFs, spliton="_")
    ORF_nams_fasta = set(list(ORF_fasta.keys()))
    gene = infiles[0][0].split("_")[0]
    for infile in infiles[1:]:
        if infile[0].startswith(gene):
            matches = infile[1]
            
    convert = pd.read_csv("%s/convert.tsv"
                          % PARAMS['sequencedir'], sep="\t", index_col=0)
    match_tab = pd.read_csv(matches, sep="\t")
    D = dict(zip(convert['id'], convert['match']))
    groups = []
    for nam in match_tab['match'].values:
        if nam in D:
            group = D[nam]
        else:
            group = "_".join(nam.split("_")[-2:])
        groups.append(group)
    match_tab['group'] = groups
    match_tab['short_id'] = match_tab['id'].str.split("_").str.get(0)
    match_tab = match_tab[match_tab['short_id'].isin(ORF_nams_fasta)]
    match_tab.to_csv(outfile, sep="\t", index=None)


@follows(mkdir("trees"))
@follows(mkdir("fastas"))
@split(makeGroups, ["fastas/*.fasta", "trees/*.tre"])
def makePhyloFastas(infiles, outfiles):
    '''
    For each of the retrovirus groups used in makeGroups, a
    fasta file is available in the phylogenies directory containing
    a set of related sequences to allow a finer classification of
    sequences in the group.
    For each group generated by makeGroups, a fasta file is built
    combining the sequences in the group with the fasta file in
    Reference_Phylogenies.
    If groups are very large (more than 40 sequences) a random
    sample of 20 sequences is used to represent the group.
    If combined groups of ERV regions and known retroviruses are very
    small, the Summary_Phylogenies sequences for the appropriate gene
    and genus are added to allow a phylogeny to be built.
    Groups are aligned using the MAFFT fftns algorithm.
    A tree is built for each group using the FastTree algorithm,
    using the -gtr and -nt options.
    An image of each tree is also generated, using the ete2 python package.
    The outputs are saved in the "fastas" directory and "trees" directory
    under the same of one known retrovirus representing the group.
    Fasta files are saved in fastas as .fasta, aligned fastas as _ali.fasta.
    Newick formatted trees are saved in the trees directory as .tre,
    png images in this directory as .png.
    '''
    for infile in infiles:
        grouptable = pd.read_csv(infile, sep="\t")
        allgroups = np.unique(grouptable['group'])
        gene = infile.split("_")[0]
        fasta = "%s_refiltered_exonerate.fasta" % gene
        PipelineERVs.makePhyloFasta(grouptable, allgroups, fasta,
                                    PARAMS['path_to_phyloseqs'],
                                    PARAMS['path_to_mafft'],
                                    PARAMS['path_to_fasttree'],
                                    PARAMS['sequencedir'])


@transform(makePhyloFastas, suffix(".tre"), ".tre")
def PhyloFastasDone(infile, outfile):
    '''
    Helper function for makePhyloFastas allowing ruffus to detect when it is
    complete.
    '''
    pass


@follows(mkdir("group_lists"))
@follows(mkdir("group_fastas"))
@follows(mkdir("summary_trees"))
@follows(mkdir("summary_fastas"))
@collate(makePhyloFastas, regex("trees/(.*?)_?([a-z]+)\_(gag|pol|env).tre"),
         r"summary_trees/\2_\3.tre")
def makeRepFastas(infiles, outfile):
    '''
    Clusters of newly identified sequences are identified in the
    trees generated by makePhyloFastas.
    For each of these clusters, a single representative sequence is selected
    These are combined with the Summary_Phylogenies sequences to build a
    single tree representing each group.  Fasta files are built containing the
    sequences, these are aligned using MAFFT fftns, trees are built with
    FastTree -gtr -nt and images are generated using the ete2 package.
    The number of sequences in the group is also shown in the phylogeny
    Each group is given a unique ID, these are shown in the phylogenies
    and the sequences in the group are listed in a .txt file in the groups
    directory.
    If a sample of sequences was taken in makePhyloFastas (for very big
    groups), the number of sequences shown in the phylogeny is corrected
    to account for this.  In these cases lists of sequences are not output
    to the groups directory.
    Fasta files of the sequences are saved in the summary_fastas directory,
    newick formatted trees and png images of trees are saved in the
    summary_trees directory.
    '''

    gene = outfile.split("_")[-1].split(".")[0]
    pp = PARAMS['path_to_phyloseqs']
    PipelineERVs.makeRepFastas(infiles, pp,
                               PARAMS['path_to_mafft'],
                               PARAMS['path_to_fasttree'],
                               PARAMS['sequencedir'],
                               outfile)


@merge(makeGroups, "summary.tsv")
def summarise(infiles, outfile):
    '''
    Statistics about the results are saved into the summary.tsv
    file in the working directory.
    A bar plot is generated showing how many ERV regions were
    identified on each chromosome for each gene, these are
    saved in the working directory as chromosome_counts_gag.png
    chromosome_counts_pol.png and chromosome_counts_env.png.
    A histogram is generated showing the distribution of maximum ORF
    lengths for each gene.  These are saved in the working directory as
    as orf_lengths_gag.png, orf_lengths_pol.png and orf_lengths_env.png.
    A pie chart is generated showing the distribution between retroviral
    genera identified in the host for each gene.  These are saved in the
    working directory as genera_gag.png, genera_pol.png and genera_env.png.
    A bar chart is generated for each gene and genus showing the
    number of ERVs assigned to each group by makeGroups.
    These are saved in the working directory as groups_gene_genus.png.
    '''
    best = infiles
    gagf = best[1]
    polf = best[2]
    envf = best[0]
    try:
        gag = pd.read_csv(gagf, sep="\t", header=0, index_col=0)
    except:
        gag = []

    try:
        pol = pd.read_csv(polf, sep="\t", header=0, index_col=0)
    except:
        pol = []

    try:
        env = pd.read_csv(envf, sep="\t", header=0, index_col=0)
    except:
        env = []

    PipelineERVs.summary(gag, pol, env, outfile,
                         PARAMS['working_directory'])


if __name__ == '__main__':
    cmdline.run(options)
