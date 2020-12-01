#!/usr/bin/env python3
def raiseError(error, *args, log=None):
    print (*args)
    err = error(*args)
    if log:
        log.error(err)
    raise (err)


def DatabaseNotFoundError(*args):
    return RuntimeError("Database fasta file %s does not exist" % args[0])


def NoIniError(*args):
    return RuntimeError("""A copy of the configuration file pipeline.ini needs to be in your working directory""")


def FastaIndexError(*args):
    return RuntimeError("""Indexing the input Fasta file was not possible, please check your input file for errors""")


def InputFileNotSpecifiedError(*args):
    return RuntimeError("Input file needs to be specified in the pipeline.ini file")


def InputFileNotFoundError(*args):
    return FileNotFoundError("Input file %s not found" % (args[0]))


def ToolNotInPathError(*args):
    if len(args) == 2:
        return RuntimeError("%s is not at the location %s specified in your pipeline.ini" % (args[0], args[1]))
    else:
        return RuntimeError("%s is not in your $PATH" % (args[0]))


def ToolPathNotSpecifiedError(*args):
    return RuntimeError("Path to %s not specified in your pipeline.ini" % args[0])


def TooManyRunsWarn(*args):
    return RuntimeWarning("""With your current settings, your input genome will be divided into %i batches and Exonerate will run %i times. Pipeline will continue as genome_force is set to True.""" % (args[0], args[1]))


def TooManyRunsError(*args):
    return RuntimeError("""With your current settings, your input genome will be divided into %i batches and Exonerate will run %i times. To continue please set genome_split to True in the pipeline.ini to reduce the number of runs or force the pipeline to run with this large number of divisions by setting genome_force to True in the pipeline.ini.""" % (args[0], args[1]))


def MissingChromsError(*args):
    return RuntimeError("""Not all chromosomes in keepchroms.txt are found in the input fasta file""")


def UBLASTDBError(*args):
    return RuntimeError("""Error generating UBLAST database for file %s - see log file""" % args[0])


def SoftwareError(*args):
    return RuntimeError("""Error running %s on %s - see log file""" % (args[0], args[1]))