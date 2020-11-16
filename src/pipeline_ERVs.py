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
import pathlib


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

    PipelineERVs.runExonerate(infiles[1], infiles[0],
                              outfile, log,
                              PARAMS['paths']['path_to_exonerate'])


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
    combined = PipelineERVs.combineBeds(beds)
    if combined is not None:
        log.info("%i records identified in combined bed file %s" % (
            len(combined), outfiles[0]))
        combined.to_csv(outfiles[0], sep="\t", index=None, header=None)

        merged = PipelineERVs.mergeBed(outfiles[0],
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
        statement = ['bedtools',
                     'getfasta',
                     '-s', '-fi', PARAMS['input']['genome'],
                     '-bed', infile]
        log.info("Generating fasta file of regions in %s: %s" % (infile,
                                                                 statement))

        P = subprocess.run(statement,
                           stdout=open(outfile, "w"),
                           stderr=subprocess.PIPE)
        if P.returncode != 0:
            log.error(P.stderr)
            err = RuntimeError("Error converting bed file to fasta - see log file")
            log.error(err)
            raise err


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
        F = PipelineERVs.makeFastaDict(infile)
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
        err = FileNotFoundError("The path %s to the usearch executable is not valid" % usearch)
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
        PipelineERVs.filterFasta(L, fasta_in, fasta_out, log, split=False)


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
        reference_ERVs = "%s/ERV_db/all_ERVS.fasta" % PARAMS[
            'database']['path_to_ERVsearch']
        PipelineERVs.classifyWithExonerate(reference_ERVs,
                                           fasta, outfiles[0],
                                           exonerate_path,
                                           exonerate_minscore,
                                           log)

    log.info("Converting raw exonerate output %s to a table" % outfiles[0])
    res = pd.read_csv(outfiles[0],
                      sep="\t", header=None, names=['id', 'match', 'score'])
    log.info("Finding the highest scoring hit for each putative ERV in \
             %s" % outfiles[0])
    res = PipelineERVs.findBestExonerate(res, gene)
    res.to_csv(outfiles[1], sep="\t", index=None)

    log.info("Generating a FASTA file for the results in %s" % outfiles[1])
    L = list(set(res['id']))
    PipelineERVs.filterFasta(L, fasta, outfiles[2], log, split=False)


@follows(mkdir("ORFs.dir"))
@transform(classifyWithExonerate,
           regex(
               "exonerate_classification.dir/(.*)_all_matches_exonerate.tsv"),
           [r"ORFs.dir/\1_orfs_raw.fasta",
            r"ORFs.dir/\1_orfs_nt.fasta",
            r"ORFs.dir/\1_orfs_aa.fasta"])
def ORFs(infiles, outfiles):
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
        PipelineERVs.runTranseq(fasta,
                                outfiles[0],
                                PARAMS['orfs']['translation_table'],
                                log)
        PipelineERVs.filterTranseq(fasta,
                                   outfiles[0], outfiles[1], outfiles[2],
                                   int(PARAMS['orfs']['min_orf_len']),
                                   PARAMS['input']['genome'], log)


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
