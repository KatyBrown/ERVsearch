from ruffus import *

import PipelineERVs as PipelineERVs
import os
import pandas as pd
import numpy as np
import shutil
import ruffus.cmdline as cmdline

PARAMS = PipelineERVs.getParameters("./pipeline.ini")
for PARAM in PARAMS:
    print "%s\t%s\n" % (PARAM, PARAMS[PARAM])

parser = cmdline.get_argparse(description='Pipeline ERVs')

options = parser.parse_args()

logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)


@follows(mkdir("UBLAST_db"))
@follows(mkdir("host_chromosomes"))
@split("%s/%s.fa" % (PARAMS['genome_directory'], PARAMS['genome']),
       r"host_chromosomes/*.fa")
def genomeToChroms(infile, outfiles):
    '''
    Splits the host genome provided by the user into one fasta file for each
    chromosome.  If the genome is not assembled into chromosomes, this is
    specified in  pipeline.ini. The pipeline.ini parameter nchroms is then
    used - the contigs or scaffolds are concatenated into nchroms fasta files.
    This input type is allows much quicker screening with Exonerate.
    This function generates a series of fasta files which are stored in the
    host_chromosomes directory.
    '''
    if int(PARAMS['has_chroms']) == 0:
        PipelineERVs.makeChroms(infile, PARAMS['n_chroms'])
    else:
        PipelineERVs.splitChroms(infile)


@transform(genomeToChroms, suffix(".fa"), ".fa.fai")
def indexChroms(infile, outfile):
    '''
    Indexes all chromosomes in the host_chromosomes directory (or chromosome
    constructs generated by genomeToChroms) using Samtools faidx,
    for fast sequence retrieval later, and saves the index files.
    '''
    fasta = infile
    temp = "%s/temp.fa" % PARAMS['working_directory']
    path_to_reformat = PARAMS['path_to_exonerate'].replace(
        "bin/exonerate", "bin/")
    os.system("%sfastareformat %s > %s" % (path_to_reformat,
                                           fasta, temp))
    shutil.move(temp, fasta)
    os.system("%s faidx %s"
              % (PARAMS['path_to_samtools'], infile))


@follows(indexChroms)
@follows(mkdir("raw_exonerate_output"))
@follows(mkdir("parsed_exonerate_output"))
@follows(genomeToChroms)
@transform("%s/*fa" % PARAMS["sequencedir"],
           regex("%s/(.*).fa" % PARAMS["sequencedir"]),
           add_inputs(genomeToChroms),
           (r"parsed_exonerate_output/\1_results.tsv",
            r"parsed_exonerate_output/\1_results.bed"))
def runandParseExonerate(infiles, outfiles):
    '''
    Runs the protein2dna algorithm in the Exonerate software package with
    the host chromosomes in host_chromosomes as target sequences and the
    ERV_Amino_Acid_DB fasta files as query sequences.
    Output is parsed to remove sequences shorter than the "overlap" parameter
    in pipeline.ini and sequences containing introns and combined into a
    tab-delimited file
    The raw output of Exonerate is stored in the raw_exonerate_output directory
    This step is carried out with low stringency as results are later filtered
    using UBLAST and Exonerate (in the runUBLASTCheck and makeGroups steps).
    The parsed output is stored in the parsed_exonerate_output directory as
    gags_results.tsv, pols_results.tsv and envs_results.tsv.
    A bed file (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    is also generated corresponding to the regions in each parsed output file -
    these are stored in parsed_exonerate_output as gags.bed, pols.bed, envs.bed
    '''
    fastaFile = infiles[0]
    chromFiles = infiles[1:]
    fastastem = fastaFile.split("/")[-1].replace(".fa", "")
    min_hit_length = int(PARAMS['min_hit_length'])
    cnames = ["query_id", "query_start", "query_end", "query_strand",
              "target_id", "target_start", "target_end",
              "target_strand", "score", "details", "length"]
    bigdf = pd.DataFrame(columns=cnames)
    for chrom in chromFiles:
        chromstem = chrom.split("/")[-1].replace(".fa", "")
        out = "raw_exonerate_output/%s_%s.tsv" % (fastastem, chromstem)

        PipelineERVs.runExonerate(fastaFile, chrom, out,
                                  PARAMS['path_to_exonerate'])

        smalldf = PipelineERVs.filterExonerate(out, min_hit_length)

        bigdf = bigdf.append(smalldf)
    bigdf.to_csv(outfiles[0], sep="\t")
    bedcols = bigdf[['target_id', 'target_start', 'target_end',
                     'target_strand', 'score', 'query_id']]
    bedcols.to_csv(outfiles[1], sep="\t", header=False, index=False)


@transform(runandParseExonerate, suffix("_results.tsv"),
           "_merged.bed")
def mergeOverlaps(infiles, outfile):
    '''
    Overlapping regions of the genome detected by Exonerate with
    similarity to the same retroviral gene are merged into
    single regions.  This is performed using bedtools merge on
    the bed files output by runAndParseExonerate.
    Merged bed files are stored in parsed_exonerate_output as
    gags_merged.bed, pols_merged.bed and envs_merged.bed.
    '''
    bedfile = infiles[1]
    os.system('%s merge -i %s > %s' % (PARAMS['path_to_bedtools'],
                                       bedfile, outfile))


@transform(mergeOverlaps, suffix("_merged.bed"),
           add_inputs(genomeToChroms),
           ".fa")
def makeFasta(infiles, outfile):
    '''
    Fasta files are generated containing the sequences
    of the merged regions of the genome identified using
    mergeOverlaps.
    These are extracted from the host chromosomes using Samtools.
    The output files are stored in parsed_exonerate_output as
    gags.fa, pols.fa and envs.fa.
    '''
    string = ""
    D = dict()
    with open(infiles[0]) as inf:
        for line in inf:
            line = line.strip().split("\t")
            chrom = line[0]
            if chrom not in D:
                D[chrom] = []
            D[chrom].append("%s:%s-%s" % (line[0], line[1], line[2]))

    for chrom in D:
        string = " ".join(D[chrom])
        os.system('%s faidx host_chromosomes/%s.fa %s >> %s'
                  % (PARAMS['path_to_samtools'], chrom, string, outfile))


@transform("%s/*fa" % PARAMS["sequencedir"],
           regex("%s/(.*).fa" % PARAMS["sequencedir"]),
           r'UBLAST_db/\1.udb')
def makeUBLASTDb(infile, outfile):
    '''
    USEARCH requires an indexed database of query sequences to run.
    This function generates this database for the three
    ERV_Amino_Acid_DB fasta files.
    '''
    os.system("%s/usearch -makeudb_ublast -quiet \
               %s -output %s --alpha aa" % (
        PARAMS['path_to_usearch'], infile, outfile))


@transform(makeFasta, suffix(".fa"), add_inputs(makeUBLASTDb),
           ("_table.tsv", "_alignments.txt", "_filtered.fa"))
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
    fastanam = fasta.split("/")[-1].replace(".fa", "")
    for inf in infiles[1:]:
        if inf.split("/")[-1].replace(".udb", "") == fastanam:
            dbnam = inf
    id_thresh = PARAMS['usearch_id']
    min_coverage = PARAMS['usearch_coverage']
    os.system("""%s/usearch -ublast %s \
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


@transform(runUBLASTCheck, suffix(".tsv"),
           add_inputs(PARAMS['path_to_refs']),
           ["_matches.tsv", "_bestmatches.tsv"])
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
    gene = infiles[0][0].split("/")[1][0:3]
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
    --score 500\
    | sort -n \
    > %s """ % (PARAMS['path_to_exonerate'],
                fasta, refs, outfiles[0]))
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
    rtab.to_csv(outfiles[1], sep="\t")


@transform(findBest, suffix("_matches.tsv"),
           add_inputs(makeFasta), ["_orfs.tsv", "_orfs.fa"])
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
    matches = infiles[0][1]
    fasta = matches.replace("_table_bestmatches.tsv", "_filtered.fa")
    PipelineERVs.getORFS(matches, fasta,
                         outfiles, PARAMS['path_to_emboss'])


@transform(ORFs, suffix("_orfs.tsv"), "_groups.tsv")
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
    matches = infiles[0]
    convert = pd.read_csv("%s/convert.tsv"
                          % PARAMS['sequencedir'], sep="\t", index_col=0)
    PipelineERVs.group(matches, convert, outfile)


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
    for filename in infiles:
        try:
            grouptable = pd.read_csv(filename, sep="\t", index_col=0)
        except:
            grouptable = []
        if len(grouptable) != 0:
            allgroups = np.unique(grouptable['group'])
            fasta = filename.replace("_table_groups.tsv", "_filtered.fa")
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
@follows(mkdir("summary_trees"))
@follows(mkdir("summary_fastas"))
@collate(PhyloFastasDone, regex("trees/(.*)\_(\S+)\_(\S+).tre"),
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
    gene = outfile.split("_")[-1].replace(".tre", "")
    wd = PARAMS['working_directory']
    genetab = pd.read_csv(
        "%s/parsed_exonerate_output/%ss_table_groups.tsv"
        % (wd, gene), sep="\t", index_col=0)

    pp = PARAMS['path_to_phyloseqs']
    PipelineERVs.makeRepFastas(infiles, genetab, wd, pp,
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
