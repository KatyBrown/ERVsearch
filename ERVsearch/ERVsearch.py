#!/usr/bin/env python3

# This is the main ruffus pipeline for ERVsearch.

from ruffus import follows, split, transform, mkdir, formatter, originate
from ruffus import collate, regex, subdivide, merge
from ruffus.combinatorics import product

import os
import pandas as pd
import numpy as np
import ruffus.cmdline as cmdline
import logging
import subprocess
import configparser
import pathlib
import math
import HelperFunctions
import Fasta
import Bed
import Exonerate
import ORFs
import Trees
import Summary
import Regions

parser = cmdline.get_argparse(description='Pipeline ERVs')
options = parser.parse_args()

# Setup to read the ini config file
PARAMS = configparser.ConfigParser()

# Check the ini config file exists
if not os.path.exists("pipeline.ini"):
    err = RuntimeError("A copy of the configuration file pipeline.ini needs\
                        to be in your working directory")
    raise (err)

# Read the parameters from the ini config file
PARAMS.read("pipeline.ini")

# Read the output file stem for the log file
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

# Setup custom gene databases
genes = []
PARAMS['gene'] = {}

# The default is to use /ERV_db/gag.fasta /ERV_db/pol.fasta and
# /ERV_db/env.fasta in the pipeline directory but the user can
# also provide their own FASTA files.

usecustom = PARAMS['database']['use_custom_db'] == "True"

for gene in ['gag', 'pol', 'env']:
    if PARAMS['database'][gene] != "None" or not usecustom:
        if usecustom:
            # read the database path from the ini file if needed
            genedb = PARAMS['database'][gene]
        else:
            # otherwise use the default database
            genedb = "%s/ERV_db/%s.fasta" % (
                PARAMS['database']['path_to_ERVsearch'], gene)
        if not os.path.exists(genedb):
            # if there's a custom database specified, check the file exists
            err = RuntimeError(
                "Database fasta file %s does not exist" % genedb)
            log.error(err)
            raise err
        # add this info to the parameter dictionary
        PARAMS['gene'][gene] = genedb
        # keep track of which genes have databases
        genes.append(gene)

# Write the parameters to file
PARAMS.write(open("%s_parameters.txt" % outstem, "w"))

# Store the plot format (as this is used in the ruffus calls)
plot_format = PARAMS['plots']['format']


@originate("init.txt")
def initiate(outfile):
    '''
    Check the environment before starting.
    Check that:
        The input file exists.
        The correct path to ERVsearch is provided.
        samtools, bedtools, FastTree and mafft are in the PATH.
        The correct paths to usearch and exonerate are provided.
    '''
    HelperFunctions.quickCheck(PARAMS, log)
    out = open(outfile, "w")
    out.write("Passed initial checks")
    out.close()


@follows(initiate)
@follows(mkdir("host_chromosomes.dir"))
@split(PARAMS['genome']['file'], r"host_chromosomes.dir/*.fasta")
def genomeToChroms(infile, outfiles):
    '''
    Splits the host genome provided by the user into one fasta file for each
    chromosome, scaffold or contig.

    An unzipped copy of zipped and gzipped fasta files will be
    created.

    This function generates a series of fasta files which are stored in the
    host_chromosomes.dir directory.
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

    statement = ["samtools", "faidx", "genome.fa"]
    log.info("Indexing fasta file: %s" % " ".join(statement))
    s = subprocess.run(statement)
    if s.returncode != 0:
        err = RuntimeError("""Indexing the input Fasta file was not possible, \
                           please check your input file for errors""")
        log.error(err)
        raise(err)

    statement = ["wc", "-l", "genome.fa.fai"]
    log.info("Counting chromosomes in fasta file genome.fa: \
             %s" % (" ".join(statement)))
    P = subprocess.run(statement, stdout=subprocess.PIPE)

    if PARAMS['genome']['split'] == "True":
        nsplits = int(PARAMS['genome']['split_n'])
        n_per_split = math.ceil(
            int(P.stdout.decode().split(" ")[0]) / nsplits)
    else:
        n_per_split = 1
    log.info("%s chromosomes will be written per file" % n_per_split)
    Fasta.splitChroms("genome.fa", log, n=n_per_split)


@follows(mkdir("gene_databases.dir"))
@originate(["gene_databases.dir/%s.fasta" % gene
            for gene in PARAMS['gene'].keys()])
def prepDBS(outfile):
    gene = os.path.basename(outfile).split(".")[0]
    genedb = PARAMS['gene'][gene]
    statement = ['cp',
                 genedb,
                 "gene_databases.dir/%s.fasta" % gene]
    log.info("Making a copy of database %s: %s" % (genedb,
                                                   " ".join(statement)))
    subprocess.run(statement)


@follows(mkdir("raw_exonerate_output.dir"))
@follows(mkdir("clean_exonerate_output.dir"))
@follows(genomeToChroms)
@product(genomeToChroms, formatter(),
         prepDBS,
         formatter(),
         r'raw_exonerate_output.dir/{basename[1][0]}_{basename[0][0]}.tsv')
def runExonerate(infiles, outfile):
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

    Exonerate.runExonerate(infiles[1], infiles[0],
                           outfile, log,
                           PARAMS['paths']['path_to_exonerate'])


@transform(runExonerate, regex("raw_exonerate_output.dir/(.*).tsv"),
           [r'clean_exonerate_output.dir/\1_unfiltered.tsv',
            r'clean_exonerate_output.dir/\1_filtered.tsv',
            r'clean_exonerate_output.dir/\1.bed'])
def cleanExonerate(infile, outfiles):
    # Minimum size number of aas in a hit t
    min_hit_length = int(PARAMS['exonerate']['min_hit_length'])
    Exonerate.filterExonerate(infile, outfiles, min_hit_length, log)


@follows(mkdir("gene_bed_files.dir"))
@collate(cleanExonerate,
         regex("clean_exonerate_output.dir/([a-z]+)_(.*)_unfiltered.tsv"),
         [r"gene_bed_files.dir/\1_all.bed",
          r"gene_bed_files.dir/\1_merged.bed"])
def mergeOverlaps(infiles, outfiles):
    '''
    Overlapping regions of the genome detected by Exonerate with
    similarity to the same retroviral gene are merged into
    single regions.  This is performed using bedtools merge on
    the bed files output by cleanExonerate.

    If there is a gap of less than PARAMS['exonerate']['overlap'] between the
    regions they will be merged.

    Merged bed files are stored as
    gags_merged.bed, pols_merged.bed and envs_merged.bed.
    '''
    log.info("Generating combined bed file %s" % outfiles[0])
    beds = [inf[2] for inf in infiles]
    combined = Bed.combineBeds(beds)
    if combined is not None:
        log.info("%i records identified in combined bed file %s" % (
            len(combined), outfiles[0]))
        combined.to_csv(outfiles[0], sep="\t", index=None, header=None)

        merged = Bed.mergeBed(outfiles[0],
                              PARAMS['exonerate']['overlap'],
                              log)
        if merged is not None:
            log.info("Writing merged bed file %s with %i lines" % (
                outfiles[1], len(merged)))
            merged.to_csv(outfiles[1], sep="\t", index=None, header=None)
        else:
            log.info("No ERVs identified in %s" % outfiles[1])
            pathlib.Path(outfiles[1]).touch()
    else:
        log.info("No records identified in combined bed file %s" % outfiles[0])
        pathlib.Path(outfiles[0]).touch()
        pathlib.Path(outfiles[1]).touch()


@follows(mkdir("gene_fasta_files.dir"))
@transform(mergeOverlaps, regex("gene_bed_files.dir/(.*)_all.bed"),
           r'gene_fasta_files.dir/\1_merged.fasta')
def makeFastas(infiles, outfile):
    '''
    Fasta files are generated containing the sequences
    of the merged regions of the genome identified using
    mergeOverlaps.
    These are extracted from the host chromosomes using beodtols.
    The output files are stored in gene_fasta_files.dir.
    '''
    infile = infiles[1]

    if os.stat(infile).st_size == 0:
        pathlib.Path(outfile).touch()
    else:
        Bed.getFasta(infile, outfile, log)


@transform(makeFastas, regex("gene_fasta_files.dir/(.*)_merged.fasta"),
           r"gene_fasta_files.dir/\1_merged_renamed.fasta")
def renameFastas(infile, outfile):
    '''
    Rename the genes in the fasta files so each record has a numbered unique
    ID (gag1, gag2 etc).
    Also removes ":" from sequence names
    '''
    if os.stat(infile).st_size == 0:
        pathlib.Path(outfile).touch()
    else:
        out = open(outfile, "w")
        F = Fasta.makeFastaDict(infile)
        log.info("Generating gene IDs for %s" % infile)
        gene = os.path.basename(infile).split("_")[0]
        for i, nam in enumerate(F):
            out.write(">%s%s_%s\n%s\n" % (gene, i+1, nam.replace(":", "-"),
                                          F[nam]))
        out.close()


@follows(mkdir("UBLAST_db.dir"))
@transform(prepDBS,
           formatter(),
           r'UBLAST_db.dir/{basename[0]}_db.udb')
def makeUBLASTDb(infile, outfile):
    '''
    USEARCH requires an indexed database of query sequences to run.
    This function generates this database for the three
    ERV_Amino_Acid_DB fasta files.
    '''
    usearch = PARAMS['paths']['path_to_usearch']
    if not os.path.exists(usearch):
        err = FileNotFoundError("The path %s to the usearch executable \
                                is not valid" % usearch)
        log.error(err)
        raise err
    statement = [usearch,
                 '-makeudb_ublast',
                 infile,
                 '-output',
                 outfile,
                 '-quiet']
    log.info("Building usearch database for %s: %s" % (infile,
                                                       " ".join(statement)))
    P = subprocess.run(statement,
                       stderr=subprocess.PIPE)
    if P.returncode != 0:
        log.error(P.stderr)
        err = RuntimeError(
            "Error making usearch database %s - see log file" % outfile)
        log.error(err)
        raise err


@follows(mkdir("ublast.dir"))
@collate((renameFastas, makeUBLASTDb), regex("(.*).dir/([a-z]+)_(.*)"),
         [r"ublast.dir/\2_UBLAST_alignments.txt",
          r"ublast.dir/\2_UBLAST.tsv",
          r"ublast.dir/\2_filtered_UBLAST.fasta"])
def runUBLASTCheck(infiles, outfiles):
    '''
    ERV regions in the fasta files generated by makeFasta
    are compared to the ERV_Amino_Acid_DB files for a second
    time, this time using USEARCH.

    This allows sequences with low similarity to known ERVs
    to be filtered out.  Similarity thresholds can be set in
    the pipeline.ini file (usearch_id, min_hit_length and usearch_coverage).

    The output filess are:
        ublast.dir/XXX_UBLAST_alignments.txt - UBLAST alignments output
        ublast.dir/XXX_UBLAST.tsv- tabular UBLAST output
        ublast.dir/XXX_BLAST.fasta - fasta file of UBLAST hits
    '''
    db, fasta_in = infiles
    alignments, tab, fasta_out = outfiles
    if os.stat(fasta_in).st_size == 0:
        pathlib.Path(alignments).touch()
        L = []
    else:

        min_id = PARAMS['usearch']['min_id'].strip()
        min_coverage = PARAMS['usearch']['min_coverage'].strip()
        min_length = PARAMS['usearch']['min_hit_length'].strip()

        statement = [PARAMS['paths']['path_to_usearch'],
                     '-ublast', fasta_in,
                     '-db', db,
                     "-query_cov", min_coverage,
                     "-id", min_id,
                     "-mincols", min_length,
                     '-top_hit_only',
                     '-evalue', "1",
                     '-blast6out', tab,
                     '-alnout', alignments,
                     '-quiet']
        log.info("Running UBLAST check on %s: %s" % (fasta_in,
                                                     " ".join(statement)))
        P = subprocess.run(statement, stderr=subprocess.PIPE)
        if P.returncode != 0:
            log.error(P.stderr)
            err = RuntimeError(
                "Error running usearch on %s - see log file" % fasta_in)
            log.error(err)
            raise err
        L = set([line.strip().split("\t")[0]
                 for line in open(tab).readlines()])
    # touch output files if nothing is found
    if len(L) == 0:
        pathlib.Path(tab).touch()
        pathlib.Path(fasta_out).touch()
    else:
        Fasta.filterFasta(L, fasta_in, fasta_out, log, split=False)


@follows(mkdir("exonerate_classification.dir"))
@transform(runUBLASTCheck, regex("ublast.dir/(.*)_UBLAST_alignments.txt"),
           [r"exonerate_classification.dir/\1_all_matches_exonerate.tsv",
            r"exonerate_classification.dir/\1_best_matches_exonerate.tsv",
            r"exonerate_classification.dir/\1_refiltered_exonerate.fasta"])
def classifyWithExonerate(infiles, outfiles):
    '''
    Runs the exonerate ungapped algorithm with each ERV region
    in the fasta files generated by makeFasta as queries and the
    all_ERVS.fasta fasta file as a target, to detect which known
    retrovirus is most similar to each newly identified ERV region.

    all_ERVS.fasta contains nucleic acid sequences for many known
    endogenous and exogenous retroviruses

    First all seqeunces are compared to the database and the raw output is
    saved as exonerate_classification.dir/XXX_all_matches_exonerate.tsv.
    Results need a score greater than PARAMS['exonerate']['min_score']
    against one of the genes of the same type (gag, pol or env) in the
    database.

    These results are then sorted and filtered to keep only the highest
    scoring hit for each putative ERV region, this is saved as
    exonerate_classification.dir/XXX_best_matches_exonerate.tsv.

    '''
    gene = os.path.basename(infiles[0]).split("_")[0]
    fasta = infiles[2]

    if os.stat(fasta).st_size == 0:
        for outfile in outfiles:
            pathlib.Path(outfile).touch()
    else:
        exonerate_path = PARAMS['paths']['path_to_exonerate']
        exonerate_minscore = PARAMS['exonerate']['min_score']
        reference_ERVs = "%s/ERV_db/all_ERVs_nt.fasta" % PARAMS[
            'database']['path_to_ERVsearch']
        Exonerate.classifyWithExonerate(reference_ERVs,
                                        fasta, outfiles[0],
                                        exonerate_path,
                                        exonerate_minscore,
                                        log)

    log.info("Converting raw exonerate output %s to a table" % outfiles[0])
    res = pd.read_csv(outfiles[0],
                      sep="\t", header=None, names=['id', 'match', 'score'])
    log.info("Finding the highest scoring hit for each putative ERV in \
             %s" % outfiles[0])
    res = Exonerate.findBestExonerate(res, gene)
    res.to_csv(outfiles[1], sep="\t", index=None)

    log.info("Generating a FASTA file for the results in %s" % outfiles[1])
    L = list(set(res['id']))
    Fasta.filterFasta(L, fasta, outfiles[2], log, split=False)


@follows(mkdir("ORFs.dir"))
@transform(classifyWithExonerate,
           regex(
               "exonerate_classification.dir/(.*)_all_matches_exonerate.tsv"),
           [r"ORFs.dir/\1_orfs_raw.fasta",
            r"ORFs.dir/\1_orfs_nt.fasta",
            r"ORFs.dir/\1_orfs_aa.fasta"])
def getORFs(infiles, outfiles):
    '''
    Finds the longest open reading frame in each of the ERV regions
    in the output table
    This analysis is performed using EMBOSS sixpack

    Raw sixpack output is saved in ORFs.dir as
    ORFs.dir/XXX_raw_orfs.fasta

    The start, end, length and sequence of each ORF are added to the
    output tables and saved in parsed_exonerate_output as
    gags_table_orfs.tsv, pols_table_orfs.tsv and envs_table_orfs.tsv.
    '''
    fasta = infiles[2]
    if os.stat(fasta).st_size == 0:
        for outfile in outfiles:
            pathlib.Path(outfile).touch()
    else:
        log.info("Looking for ORFs in %s" % fasta)
        ORFs.runTranseq(fasta,
                        outfiles[0],
                        PARAMS['orfs']['translation_table'],
                        log)
        ORFs.filterTranseq(fasta,
                           outfiles[0], outfiles[1], outfiles[2],
                           int(PARAMS['orfs']['min_orf_len']),
                           "genome.fa",
                           PARAMS['orfs']['translation_table'], log)


@follows(mkdir("ublast_orfs.dir"))
@collate((getORFs, makeUBLASTDb),
         regex("(.*).dir/([a-z]+)_(.*)"),
         [r"ublast_orfs.dir/\2_UBLAST_alignments.txt",
          r"ublast_orfs.dir/\2_UBLAST.tsv",
          r"ublast_orfs.dir/\2_filtered_UBLAST_nt.fasta",
          r"ublast_orfs.dir/\2_filtered_UBLAST_aa.fasta"])
def checkORFsUBLAST(infiles, outfiles):
    '''
    ERV ORFs in the fasta files generated by the ORFs function
    are compared to the ERV_Amino_Acid_DB files using USEARCH.

    This allows sequences with low similarity to known ERVs
    to be filtered out.  Similarity thresholds can be set in
    the pipeline.ini file (usearch min_id, min_hit_length and min_coverage).

    The output filess are:
        ublast_orfs.dir/XXX_UBLAST_alignments.txt - raw UBLAST alignments\
        output
        ublast_orfs.dir/XXX_UBLAST.tsv- tabular UBLAST output
        ublast_orfs.dir/XXX_BLAST.fasta - fasta file of UBLAST hits
    '''
    db = infiles[0]
    fasta_nt, fasta_aa = infiles[1][1:]
    alignments, tab, fasta_out_nt, fasta_out_aa = outfiles
    if os.stat(fasta_nt).st_size == 0:
        pathlib.Path(alignments).touch()
        L = []
    else:

        min_id = PARAMS['usearch']['min_id'].strip()
        min_coverage = PARAMS['usearch']['min_coverage'].strip()
        min_length = PARAMS['usearch']['min_hit_length'].strip()

        statement = [PARAMS['paths']['path_to_usearch'],
                     '-ublast', fasta_nt,
                     '-db', db,
                     "-query_cov", min_coverage,
                     "-id", min_id,
                     "-mincols", min_length,
                     '-top_hit_only',
                     '-evalue', "1",
                     '-blast6out', tab,
                     '-alnout', alignments,
                     '-quiet']
        log.info("Running UBLAST check on %s: %s" % (fasta_nt,
                                                     " ".join(statement)))
        P = subprocess.run(statement, stderr=subprocess.PIPE)
        if P.returncode != 0:
            log.error(P.stderr)
            err = RuntimeError(
                "Error running usearch on %s - see log file" % fasta_nt)
            log.error(err)
            raise err
        L = set([line.strip().split("\t")[0]
                 for line in open(tab).readlines()])
    # touch output files if nothing is found
    if len(L) == 0:
        pathlib.Path(tab).touch()
        pathlib.Path(fasta_out_nt).touch()
        pathlib.Path(fasta_out_aa).touch()
    else:
        Fasta.filterFasta(L, fasta_nt, fasta_out_nt, log, split=False)
        Fasta.filterFasta(L, fasta_aa, fasta_out_aa, log, split=False)


@follows(mkdir("grouped.dir"))
@transform(checkORFsUBLAST,
           regex("ublast_orfs.dir/(.*)_UBLAST_alignments.txt"),
           r"grouped.dir/\1_groups.tsv")
def assignGroups(infiles, outfile):
    '''
    The retroviruses in all_ERVs_Fasta have been classified
    into groups based on sequence similarity.

    Each group is named after a single representative ERV.

    The newly identified ERV regions are classified into the
    same groups based on the most similar sequence identified with
    UBLAST.

    The assigned group is added to the output tables, these are
    saved as gags_table_groups.tsv, pols_table_groups.tsv and
    envs_table_groups.tsv.
    '''
    ORF_file = infiles[2]
    ORF_fasta = Fasta.makeFastaDict(ORF_file, spliton="_")

    # this table has the group information - read it and put it into
    # a dictionary
    convert = pd.read_csv("%s/ERV_db/convert.tsv"
                          % PARAMS['database']['path_to_ERVsearch'],
                          sep="\t")
    D = dict(zip(convert['id'], convert['match']))

    log.info("Assigning %s to a group" % infiles[2])
    # Read the UBLAST matches - these should already be the best UBLAST
    # match
    match_tab = pd.read_csv(infiles[1], sep="\t", header=None)
    match_tab = match_tab[[0, 1, 2, 3, 10, 11]]
    match_tab.columns = ['name', 'match', 'perc_identity', 'alignment_length',
                         'evalue', 'bit_score']
    match_tab['ID'] = [x.split("_")[0] for x in match_tab['name']]

    segs = [ORFs.splitNam(nam) for nam in match_tab['name']]
    stab = pd.DataFrame(segs)

    match_tab = match_tab.merge(stab)

    groups = []

    for nam in match_tab['match'].values:
        nam_stem = "_".join(nam.split("_")[:-2])
        if nam_stem in D:
            group = D[nam_stem]
        else:
            group = "_".join(nam_stem.split("_")[-2:])
        groups.append(group)
    match_tab['group'] = groups
    match_tab['evalue'] = ["%.3e" % e for e in match_tab['evalue']]
    match_tab['genus'] = match_tab['match'].str.split("_").str.get(-4)
    match_tab = match_tab[match_tab['ID'].isin(ORF_fasta)]
    match_tab.to_csv(outfile, sep="\t", index=None)


@follows(assignGroups)
def Screen():
    pass


@follows(mkdir("summary_tables.dir"))
@follows(mkdir("summary_plots.dir"))
@merge(assignGroups,
       [r"summary_tables.dir/exonerate_initial_summary.txt",
        r"summary_plots.dir/exonerate_initial_lengths.%s" % plot_format,
        r"summary_plots.dir/exonerate_initial_by_sequence.%s" % plot_format])
def summariseScreen(infiles, outfiles):
    Summary.summariseExonerateInitial(infiles, outfiles, log, genes,
                                      PARAMS['plots'])
    Summary.summariseUBLAST(infiles, outfiles, log, genes,
                            PARAMS['plots'])


@follows(mkdir("group_fastas.dir"))
@subdivide(assignGroups, regex("grouped.dir/([a-z]+)_groups.tsv"),
           [r"group_fastas.dir/\1_*.fasta",
            r"group_fastas.dir/\1_*_A.fasta"])
def makeGroupFastas(infile, outfiles):
    '''
    Two sets of reference fasta files are available (files are stored in
        phylogenies/group_phylogenies and phylogenies/summary_phylogenies)

        group_phylogenies - groups of closely related ERVs for fine
        classification of sequences
        summary_phylogenies - groups of most distant ERVs for broad
        classification of sequences

    For each of the groups assigned in makeGroups, if there are more than 8
    members of the group a fasta file will
    be created with the other members of the same group from group_phylogenies.
    Otherwise, a fasta file will be created with the summary_phylogenies
    group for this gene and genus.
    '''
    grouptable = pd.read_csv(infile, sep="\t")
    gene = os.path.basename(infile).split("_")[0]
    fasta = "ublast_orfs.dir/%s_filtered_UBLAST_nt.fasta" % gene
    Trees.makeGroupFasta(grouptable, fasta,
                         PARAMS['database']['path_to_ERVsearch'],
                         log)


@follows(mkdir("group_trees.dir"))
@transform(makeGroupFastas, regex("group_fastas.dir/(.*)_A.fasta"),
           r"group_trees.dir/\1.tre")
def makeGroupTrees(infile, outfile):
    '''
    A tree is built for each group using the FastTree algorithm,
    using the -gtr and -nt options.
    '''
    Trees.buildTree(infile, outfile, log)


@transform(makeGroupTrees, regex("group_trees.dir/(.*).tre"),
           r"group_trees.dir/\1.%s" % PARAMS['trees']['format'])
def drawGroupTrees(infile, outfile):
    '''
    An image of each tree is  generated, using the ete3 python package.
    '''
    if PARAMS['trees']['use_gene_colour'] == "True":
        hlcolour = PARAMS['plots']['%s_colour' % (os.path.basename(
            infile).split("_")[0])]
    else:
        hlcolour = PARAMS['trees']['highlightcolour']
    Trees.drawTree(infile, outfile, PARAMS['trees']['maincolour'],
                   hlcolour,
                   PARAMS['trees']['outgroupcolour'],
                   PARAMS['trees']['dpi'])


@follows(mkdir("summary_fastas.dir"))
@follows(mkdir("group_lists.dir"))
@collate((makeGroupFastas, makeGroupTrees),
         regex("group_[a-z]+.dir/([a-z]+)(_?.*?)_([a-z]+)\.([a-z]+)"),
         [r"summary_fastas.dir/\1_\3.fasta",
          r"summary_fastas.dir/\1_\3_A.fasta"])
def makeSummaryFastas(infiles, outfiles):
    '''
    Clusters of newly identified sequences are identified in the
    trees generated by makePhyloTrees

    For each of these clusters, a single representative sequence is selected

    These are combined with the group_phylogenies and summary_phylogenies
    sequences to build a single tree representing each group.

    Fasta files of the sequences are saved in the summary_fastas.dir directory,

    '''
    inf = np.array(infiles)
    fastas = np.sort(inf[np.char.endswith(inf, ".fasta")])
    trees = np.sort(inf[np.char.endswith(inf, ".tre")])
    Trees.makeRepFastas(fastas, trees,
                        PARAMS['database']['path_to_ERVsearch'],
                        outfiles, log)


@follows(mkdir("summary_trees.dir"))
@transform(makeSummaryFastas,
           regex("summary_fastas.dir/(.*).fasta"),
           r"summary_trees.dir/\1.tre")
def makeSummaryTrees(infiles, outfile):
    Trees.buildTree(infiles[1], outfile, log)


@transform(makeSummaryTrees,
           regex("summary_trees.dir/(.*).tre"),
           r"summary_trees.dir/\1.%s" % PARAMS['trees']['format'])
def drawSummaryTrees(infile, outfile):
    if PARAMS['trees']['use_gene_colour'] == "True":
        hlcolour = PARAMS['plots']['%s_colour' % (os.path.basename(
            infile).split("_")[0])]
    else:
        hlcolour = PARAMS['trees']['highlightcolour']
    Trees.drawTree(infile, outfile, PARAMS['trees']['maincolour'],
                   hlcolour,
                   PARAMS['trees']['outgroupcolour'],
                   PARAMS['trees']['dpi'], sizenodes=True)


@follows(drawGroupTrees)
@follows(drawSummaryTrees)
def Groups():
    pass


@follows(Groups)
@merge(drawSummaryTrees,
       [r"summary_tables.dir/trees_summary.txt"])
def summariseGroups(infiles, outfiles):
    pass


@follows(mkdir("clean_beds.dir"))
@transform(assignGroups, regex("grouped.dir/(.*)_groups.tsv"),
           r"clean_beds.dir/\1.bed")
def makeCleanBeds(infile, outfile):
    df = pd.read_csv(infile, sep="\t")
    cols = HelperFunctions.getBedColumns()
    df['start'][df['strand'] == "+"] -= 1
    df = df[cols]
    df.to_csv(outfile, sep="\t", header=None, index=None)


@follows(mkdir("clean_fastas.dir"))
@transform(makeCleanBeds, regex("clean_beds.dir/(.*).bed"),
           r"clean_fastas.dir/\1.fasta")
def makeCleanFastas(infile, outfile):
    if os.stat(infile).st_size == 0:
        pathlib.Path(outfile).touch()
    else:
        Bed.getFasta(infile, outfile, log)


@follows(mkdir("ERV_regions.dir"))
@merge(makeCleanBeds,
       ["ERV_regions.dir/all_ORFs.bed",
        "ERV_regions.dir/all_regions.bed",
        "ERV_regions.dir/multi_gene_regions.bed",
        "ERV_regions.dir/regions.fasta"])
def findERVRegions(infiles, outfiles):
    combined = Bed.combineBeds(infiles)
    combined.to_csv(outfiles[0], sep="\t", index=None, header=None)
    merged = Bed.mergeBed(outfiles[0], overlap=PARAMS['regions']['maxdist'],
                          log=log, mergenames=False)
    merged.to_csv(outfiles[1], sep="\t", index=None, header=None)
    merged = merged[merged[3].str.find(",") != -1]
    merged.to_csv(outfiles[2], sep="\t", index=None, header=None)
    Bed.getFasta(outfiles[2], outfiles[3], log)


@merge(findERVRegions,
       ["ERV_regions.dir/ERV_regions_final.tsv",
        "ERV_regions.dir/ERV_regions_final.bed",
        "ERV_regions.dir/ERV_regions_final.fasta"])
def makeRegionTables(infiles, outfiles):
    results = Regions.getRegions(infiles[2], genes)
    results.to_csv(outfiles[0], sep="\t", index=None)
    cols = HelperFunctions.getBedColumns()
    cols = cols[:4] + cols[5:]
    results[cols].to_csv(outfiles[1], sep="\t", index=None, header=None)
    Bed.getFasta(outfiles[1], outfiles[2], log)


@follows(makeRegionTables)
def ERVRegions():
    pass


@follows(ERVRegions)
@merge(makeRegionTables,
       [r"summary_tables.dir/erv)regions_summary.txt"])
def summariseERVRegions(infiles, outfiles):
    pass


@follows(summariseScreen)
@follows(summariseGroups)
@follows(summariseERVRegions)
def Summarise():
    pass


@follows(Screen)
@follows(Trees)
@follows(ERVRegions)
@follows(Summarise)
def full():
    pass


if __name__ == '__main__':
    cmdline.run(options)
