#!/usr/bin/env python
# Program : qnrpipeline.py
# Author  : Fredrik Boulund
# Creation date: 2010-07-26
# Release date: 2012-07-20
# Dependencies:
#   blastclust
#   fluff.py(c)
#   HMMER
#   MAFFT
#   cdbfasta / cdbyank
#
# This file is the main component of a pipeline that searches databases using
# HMMER to find sequences matching a hidden Markov model. These sequences are
# then clustered and aligned against each other to simplify identification of
# new gene variants and maybe even new genes entirely.
#
# Copyright (C) 2012 Fredrik Boulund
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

from sys import argv, exit
from os import path, makedirs, system
from datetime import date
from math import ceil, floor
import time
import pickle
import shlex, subprocess
from optparse import OptionParser, OptionGroup

import fluff
from hmmsearch import HMMSearch
from parser import Parser
from logger import Logger
from blastclust import BLASTClusterer

ver = "QNR-search pipeline, version 0.8067 BETA" # 2012-07-20
fill_length = int(floor((78-len(ver))/2))
desc = "-"*fill_length+ver+"-"*fill_length+"""
Input database(s) in FASTA format and the pipeline will search them with HMMer's "hmmsearch" using the model specified.
This pipeline is dependent on the following external parts:
HMMER, BLASTclust, MAFFT, fluff.py(c), cdbfasta/cdbyank.
------------------------------------------------------------------------------
If you use this pipeline in your research, please cite:
Boulund et al. 2012.
The project website is located at http://bioinformatics.math.chalmers.se/qnr/
------------------------------------------------------------------------------
"""

parser = OptionParser(usage=" %prog [options] DATABASE(s)\n"\
                            "examples:\n"\
                            " %prog -M ~/model.hmm ~/database1.fasta ~/database2.fasta \n"\
                            " %prog -EL ~/hmmsearchresults/database1.fasta.hmmsearched--2010-07-04 -C 0.25 -P 85",
                            description=desc, version=ver)

# Create the standard options and add the to the OptionParser
parser.add_option("-n", "--numcpu", dest="numcpu",metavar="N",
                 help="integer - the number of CPUs, N, to use for programs with that ability [default: N=%default]")
parser.add_option("-M", "--hmm", dest="model",
                 help="the MODEL to use in hmmsearch [default: %default]")
parser.add_option("-H", "--hmmsearchonly", action="store_true", dest="hmmsearch",
                 help="boolean - only run hmmsearch on the databases with the model specified and save results to file. Can be combined with -E, -L or -A (e.g. -HEL). [default: %default]")
parser.add_option("-E", "--extracthitsonly", action="store_true", dest="extracthits",
                 help="boolean - only extract sequences from hmmsearch output that classify as potential hits and then quit, writing results to disk for inspection. Can be combined with -H, -L and -A (e.g. -EL). [default: %default]")
parser.add_option("-L", "--clusteringonly", action="store_true", dest="blastclust",
                 help="boolean - only perform cLustering of previously stored results, requires retrieved_sequences.fasta and pickled.hsseq. Can be combined with option -E or -A (e.g. -ELA). [default: %default]")
parser.add_option("-A", "--alignmentonly", action="store_true", dest="alignment",
                 help="boolean - only perform alignment of sequences within previously created clusters. Can be combined with option -L (e.g. -LA). [default: %default]")
parser.add_option("-p", "--percentidentity", dest="percent_identity",metavar="PI",
                 help="integer - the percent identity, PI, with which blastclust will perform the clustering, range 3-100 [default: PI=%default]")
parser.add_option("-C", "--coveragethreshold", dest="cov_threshold",metavar="CT",
                 help="float - the length coverage threshold, CT, used in blastclust, range 0.1-0.99 [default: CT=%default]")
parser.add_option("-k", "--classifyK", dest="classifyK", metavar="k", type="float",
                 help="float - modify the classification function parameter slope [default: k=%default]")
parser.add_option("-m", "--classifyM", dest="classifyM", metavar="m", type="float",
                 help="float - modify the classification function parameter intercept [default: m=%default]")
parser.add_option("-c", "--classifyC", dest="classifyC", metavar="c", type="float",
                 help="float - modify the classification function parameter long sequence fixed cutoff [default: c=%default]")
parser.add_option("-d", "--classifyD", dest="classifyD", metavar="d", type="float",
                 help="float - modify the classification function parameter definition of long sequence (i.e. after what fragment length to use fixed cutoff) [default: d=%default]")
parser.add_option("-a", "--addrefseq", dest="addrefseq", metavar="PATH",
                 help="add reference database with sequences in FASTA format at PATH. These sequences are added before clustering to improve/simplify cluster analysis. [default: not used]")

# Advanced options group
advanced = OptionGroup(parser, "Advanced options",
                    "Caution: some of the options might negatively affect pipeline performance!")
advanced.add_option("-s", "--noheuristics", action="store_true", dest="noheuristics",
                 help="boolean - turn hmmsearch HEURISTICS OFF (max Sensitivity) (not recommended) [default: %default]")
advanced.add_option("-D", "--retrdb", action="store_true", dest="retrdb",
                 help="boolean - retrieve sequences from original database and not from hmmsearch output file. Slower and not 'safe' for clustering or alignment. Requires a cdbfasta index file for each database [default: %default]")
advanced.add_option("-l", "--extendleft", dest="extendleft",
                 help="integer - number of residues/nucleotides to extend the HMMER aligned domain with to the left. Requires a cdbfasta index file for each database [default: not used]")
advanced.add_option("-r", "--extendright", dest="extendright",
                 help="integer - number of residues/nucleotides to extend the HMMER aligned domain with to the right. Requires a cdbfasta index file for each database [default: not used]")
advanced.add_option("--minlength", dest="minlength",
                 help="integer - minimum fragment length allowed, anything below this will be rejected by the classifier [default: %default]")
advanced.add_option("--minscore", dest="minscore",
                 help="integer - minimum fragment bit score from HMMER3, anything below this will be excluded when parsing the HMMER3 output [default: %default]")

# Developer options group
developer = OptionGroup(parser, "Developer options, use at own risk!",
                    "Caution: some of the options might negatively affect program execution and/or are generally not properly tested!")
developer.add_option("--hmmsearch_outdir", dest="hmmsearch_outdir", metavar="OUTDIR",
                 help="modify OUTput DIRectory for hmmsearch [default: %default]")
developer.add_option("--resdir", dest="resdir",
                 help="modify the DIRectory where RESulting cluster files are written [default: %default]")
developer.add_option("--alignseqpath", dest="alignseqpath",
                 help="path to fasta file with sequences to align the clusters against. It is recommended that there are not more than a couple of sequences in this file. The subroutine that performs the parsing is hardcoded to take the five PMQR-variants! [default: %default]")

# Add the advanced and developer options to the OptionParser
parser.add_option_group(advanced)
parser.add_option_group(developer)

# Set default values to all options that require them
parser.set_defaults(numcpu = 4,
                    model = "model.hmm",
                    hmmsearch_outdir = "./hmmsearchresults/",
                    noheuristics = False,
                    resdir = "./results_clusters/",
                    hmmsearch = False,  # "inverted"
                    extracthits = False, # "inverted"
                    blastclust = False, # "inverted"
                    alignment = False, # "inverted"
                    percent_identity = "90",
                    cov_threshold = "0.25",
                    alignseqpath = False,
                    retrdb = False,
                    addrefseq = False,
                    classifyK = 0.7778,
                    classifyC = 109.64,
                    classifyM = -7.89,
                    classifyD = 150.64,
                    extendleft = 0,
                    extendright = 0,
                    minlength = 25,
                    minscore = 0)

options, args = parser.parse_args() # Parse the arguments


# Check if there were any command line arguments and exit if there weren't
if len(argv)<2:
    parser.print_help()
    exit()


# Create a logfile to which all messages are written
t = time.asctime(time.localtime())

logfile = Logger('qnrsearch.log')
logfile.open()

print "\n     -----=====[ ",parser.version," ]=====-----\n\nAssigned work:"
logfile.write("\n     -----=====[  "+parser.version+"  ]=====-----\n\nAssigned work:\n")


## INITIALIZATION, PATH CREATION ETC - - - - - - - - - - - - - - - - - - - - - -
# This is a very big messy pile of if, elif, else statements... As long as it
# works as intended, I don't feel like refactoring it :).
# Only run HMMsearch. -H
if options.hmmsearch and not(options.extracthits or options.blastclust or options.alignment):
    action = "Only perform hmmsearch"

# BOTH run HMMsearch and extract hits, don't cluster or align. -HE
elif (options.hmmsearch and options.extracthits) and not(options.blastclust or options.alignment):
    action = "Perform hmmsearch and extract hits"

# Extract hits from hmmsearch AND run clustering, do not align clusters. -EL
elif (options.extracthits and options.blastclust) and not(options.hmmsearch or options.alignment):
    action = "Extract hits from hmmsearch output AND cluster the sequences AND align clusters"

# HMMsearch, extract hits and cluster, don't align. -HEL
elif (options.hmmsearch and options.extracthits and options.blastclust) and not(options.alignment):
    action = "Perform hmmsearch, extract hits and cluster"

# Extract hits from hmmsearch AND run clustering AND align clusters. -ELA
elif (options.extracthits and options.blastclust and options.alignment) and not(options.hmmsearch):
    action = "Extract hits from hmmsearch output AND cluster the sequences"

# Only extract hits from hmmsearch output. -E
elif options.extracthits and not(options.hmmsearch or options.blastclust or options.alignment):
    action = "Only extract hits from hmmsearch output"

# Only run clustering on pre-existing files. -L
elif options.blastclust and not(options.hmmsearch or options.extracthits or options.alignment):
    action = "Only performing clustering"

# Only run clustering on pre-existing files AND align clusters. -LA
elif (options.blastclust and options.alignment) and not(options.hmmsearch or options.extracthits):
    action = "Only performing clustering and cluster alignment"

# Only run alignment of pre-existing clusters. -A
elif options.alignment and not(options.hmmsearch or options.extracthits or options.blastclust):
    action = "Only performing cluster alignment"

# If no "only" options were set, run entire pipeline! -HELA
elif (not(options.hmmsearch) and not(options.blastclust) and not(options.extracthits) and not(options.alignment)) or (options.hmmsearch and options.blastclust and options.extracthits and options.alignment):
    action = "Run hmmsearch, parse its output, cluster the results and align clusters (entire pipeline)"
    options.blastclust = True
    options.extracthits = True
    options.alignment = True
    options.hmmsearch = True

print action + "\n"
logfile.write(action + "\n")

if options.hmmsearch or options.extracthits:
    if options.hmmsearch:
        message = "ERROR! Needs filename(s) for database(s) to search"
    else:
        message = "ERROR! Needs filename(s) for hmmsearch output file(s) to search"
    # Check for filename arguments
    if args == ['-'] or len(args) == 0:
        logfile.write(message + "\n")
        parser.print_help()
        exit(1)

if options.hmmsearch:
    if not path.isdir(path.abspath(options.hmmsearch_outdir)):
        makedirs(options.hmmsearch_outdir)
        logfile.write(options.hmmsearch_outdir+" created...\n")
    else:
        logfile.write(options.hmmsearch_outdir+" already exists, possibly overwriting contents, continuing...\n")

if options.alignment or not options.hmmsearch:
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        logfile.write(options.resdir+" created...\n")
    else:
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")



##---------------------------------------------------------------------------##
## - - - - - - - - - - - - CONSTANTS ETC - - - - - - - - - - - - - - - - - - ##
##---------------------------------------------------------------------------##

# Create a folder for temporary storage of lots of files.
# The idea is to clean up the uninteresting temporary files in this
# directory after usage and maybe get rid of it all together.
# It remains here for convenient debuggning and access to pipeline
# files for manual inspection/further work.
TMPDIR = "./pipeline_data/"
if not path.isdir(path.abspath(TMPDIR)):
    makedirs(TMPDIR)
    logfile.write(TMPDIR+" created...\n")
else:
    logfile.write(TMPDIR+" already exists, possibly overwriting contents, continuing...\n")

# Make sure that the path to the reference sequences is a complete path
if options.addrefseq:
    if options.addrefseq.startswith("~"):
        QNR_REFERENCE_SEQUENCES_PATH = path.expanduser(options.addrefseq)
    else:
        QNR_REFERENCE_SEQUENCES_PATH = path.abspath(options.addrefseq)


RETR_SEQ_FILEPATH = path.abspath(''.join([TMPDIR,"retrieved_sequences.fasta"]))
RESDIR = path.relpath(options.resdir)


##---------------------------------------------------------------------------##
##                  GENE FRAGMENT CLASSIFICATION FUNCTION                    ##
##---------------------------------------------------------------------------##
classificationfunction = lambda L: options.classifyK*L + options.classifyM


##---------------------------------------------------------------------------##
##                   PRINT SETTINGS TO LOGFILE AND STDOUT                    ##
##---------------------------------------------------------------------------##
settings_string = ''.join(["\nThe pipeline was started with the following settings:\n",
                           "numcpu           : ", str(options.numcpu),"\n",
                           "hmmsearch        : ", str(options.hmmsearch),"\n"
                           "extracthits      : ", str(options.extracthits),"\n",
                           "blastclust       : ", str(options.blastclust),"\n",
                           "alignment        : ", str(options.alignment),"\n",
                           "model            : ", str(options.model),"\n",
                           "heuristics       : ", str(not(options.noheuristics)),"\n",
                           "hmmsearch_outdir : ", str(options.hmmsearch_outdir),"\n",
                           "resdir           : ", str(options.resdir),"\n",
                           "percent_identity : ", str(options.percent_identity),"\n",
                           "cov_threshold    : ", str(options.cov_threshold),"\n",
                           "retrdb           : ", str(options.retrdb),"\n",
                           "addrefseq        : ", str(options.addrefseq),"\n",
                           "align_clusters   : ", str(options.alignment),"\n",
                           "alignseqpath     : ", str(options.alignseqpath),"\n",
                           "classifyC        : ", str(options.classifyC),"\n",
                           "classifyK        : ", str(options.classifyK),"\n",
                           "classifyM        : ", str(options.classifyM),"\n",
                           "classifyD        : ", str(options.classifyD),"\n",
                           "extendleft       : ", str(options.extendleft),"\n",
                           "extendright      : ", str(options.extendright),"\n",
                           "minlength        : ", str(options.minlength),"\n"])
logfile.write(settings_string)

if len(args) > 1:
    printstring = "\nThe following database files were entered at command line:\n"+'\n'.join(args)
else:
    printstring = "\nThe following database file was entered at command line:\n"+'\n'.join(args)

logfile.write(printstring+"\n")
logfile.flush()



logfile.line()



#                       ____________                                           #
#                      #            #                                          #
#                      #            #                                          #
#                       |          |                                           #
#                       |          |                                           #
#                       |          |                                           #
#                       |          |                                           #
## -------------------------------------------------------------------------- ##
## BEGIN PIPELINE - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##
##                                                                            ##

    ##---------------------------------------------------------------------------##
    ## Run hmmsearch on the entered files, assuming they are databases in fasta  ##
    ## format.                                                                   ##
    ##---------------------------------------------------------------------------##

if options.hmmsearch:
    hmms = HMMSearch(logfile)
    args = hmms.search(
                path.abspath(options.model), # Retrieve the path to the model from user set variables above
                options.hmmsearch_outdir,
                options.numcpu,
                True,
                options.noheuristics,
                args)
    # Note the change in usage of variable 'args'! It now contains the hmmsearch outputfile paths
else:
    logfile.write("Not running hmmsearch\n")
    logfile.line()


## ---------------------------------------------------------------------------- ##
#                       |          |         ____________                        #
#                       |          |        #            #                       #
#                       |          |        #            #                       #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
## ---------------------------------------------------------------------------- ##



if options.extracthits:
    ##---------------------------------------------------------------------------##
    ## Parse each of the hmmsearch output files to retrieve sequences with domain##
    ## score above MIN_SCORE that classify as interesting according to the       ##
    ## classification function                                                   ##
    ##---------------------------------------------------------------------------##

    # Note that 'args' might have changed contents from the input arguments,
    # depending of if the entire two first parts of the pipeline were run together.
    parser = Parser(logfile, classificationfunction)
    parser.parse_files(args, RETR_SEQ_FILEPATH, options.minscore, options.retrdb, options.classifyC, options.classifyD, options.extendleft, options.extendright)
else:
    logfile.write("Not extracting hmmsearch hits\n")
    logfile.line()


## ---------------------------------------------------------------------------- ##
#                       |          |         ____________                        #
#                       |          |        #            #                       #
#                       |          |        #            #                       #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
## ---------------------------------------------------------------------------- ##


if options.blastclust:

    aligner = BLASTClusterer(logfile)

    ##---------------------------------------------------------------------------##
    ##     Convert sequences to a BLAST database and cluster using blastclust    ##
    ##---------------------------------------------------------------------------##
    # If reference sequences are to be added before clustering, add them
    if options.addrefseq:
        if path.isfile(options.addrefseq):
            append_refseq = "cat "+options.addrefseq+" >> "+RETR_SEQ_FILEPATH
            system(append_refseq)
            logfile.write("Added reference sequences from "+options.addrefseq+" to set to cluster\n")

    clusters, parsedblastclust, scores_ids = aligner.run(RETR_SEQ_FILEPATH, options.numcpu, options.percent_identity, options.cov_threshold)


    # Output the identified cluster to files,
    # one with clean clusters and one with scores
    clusterfilename = RESDIR+"/identified_clusters"
    withscoresfilename = RESDIR+"/identified_clusters.scores"
    logfile.write("The identified clusters are written to: "+clusterfilename+"\n")
    clusterout = open(clusterfilename,"w")
    withscores = open(withscoresfilename,"w")
    for cluster in clusters:
        for seqID in cluster:
            clusterout.write(''.join([seqID," "]))
            for database in scores_ids:
                for info in database:
                    seqscore, domscore, seqid = info
                    if seqid in seqID:
                        withscores.write(''.join([seqID,"--",
                                         str(seqscore),"--",str(domscore),' ']))
        clusterout.write("\n")
        withscores.write("\n")
    clusterout.close()
    withscores.close()

    t = time.asctime(time.localtime())
    logfile.write("Finished clustering sequences at: "+t+"\n")
    logfile.line()
    logfile.flush()

    # pickle the clusters and scores_ids for others scripts to use as a convenience
    data = (parsedblastclust,scores_ids)
    pickle_filename = ''.join([TMPDIR,"pickled.clusters"])
    outpickle = open(pickle_filename,'wb')
    pickle.dump(data,outpickle)
    outpickle.close()
else:
    logfile.write("Not clustering sequences\n")
    logfile.line()



## ---------------------------------------------------------------------------- ##
#                       |          |         ____________                        #
#                       |          |        #            #                       #
#                       |          |        #            #                       #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
#                       |          |         |          |                        #
## ---------------------------------------------------------------------------- ##


if options.alignment:
    ##---------------------------------------------------------------------------##
    ## PERFORM MULTIPLE ALIGNMENT WITHIN CLUSTERS                                ##
    ##---------------------------------------------------------------------------##
    # This function writes to disk on its own initative so no further action
    # is required here. The function runs MAFFT on all clusters with more
    # than one member and outputs all alignments to file.
    # Additionally it produces alignments against the reference qnr-genes
    # with all clusters. Files are outputted into the RESDIR directory that
    # is assumed to exist. Note that the function contains no safeguard against
    # a nonexisting directory, but will probably crash with an error message.

    # Unpickle parsedblastclust (needed for being able to run alignment separately)
    # Really unnecessary to do every time but does not really matter
    try:
        parsedblastclust, scores_ids = pickle.load(open(''.join([TMPDIR,"pickled.clusters"]),'rb'))
    except IOError:
        logfile.write("Could not read the pickled clusters (pickled.clusters)\n")
        fluff.cleanup(TMPDIR)
    except pickle.UnpicklingError:
        logfile.write("Could not unpickle pickled.clusters!\n")
        fluff.cleanup(TMPDIR)

    try:
        # If QNR_REFERENCE_SEQUENCES_PATH is invalid, do not align against
        # reference sequences and notify user path was invalid.
        if options.alignseqpath:
            # Log that we started aligning
            t = time.asctime(time.localtime())
            logfile.write(t+"\nPerforming multiple alignment within identified clusters and "\
                  "cluster vs reference genes\n")
            logfile.flush()

            # Check that the reference sequence path is valid
            if path.isfile(QNR_REFERENCE_SEQUENCES_PATH):
                retcode = fluff.malign_clusters(parsedblastclust, RESDIR,
                                                QNR_REFERENCE_SEQUENCES_PATH,
                                                TMPDIR)
            else:
                logfile.write("Filename supplied with -R flag was invalid, NOT aligning"
                       " against any reference sequences.\n")
                retcode = fluff.malign_clusters(parsedblastclust, RESDIR, "", TMPDIR)
            if retcode == 1:
                logfile.write("No clusters contained more than one member, nothing to align!\n")
            elif retcode == 2:
                logfile.write("Could not align sequences using 'mafft'. Is MAFFT properly installed?\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            else:
                logfile.write("Multiple alignments complete\n" \
                       "MAFFT alignments are in files " \
                       "'*.aligned' in "+RESDIR+"\n")
        else:
            # Log that we started aligning
            t = time.asctime(time.localtime())
            logfile.write(t+"\nPerforming multiple alignment within identified clusters\n")
            logfile.flush()

            # Align the clusters
            retcode = fluff.malign_clusters(parsedblastclust, RESDIR,
                                            refseqpath="",
                                            workdir=TMPDIR)
            if retcode == 1:
                logfile.write("No clusters contained more than one member, nothing to align!\n")
            elif retcode == 2:
                logfile.write("Could not align sequences using 'mafft'. Is MAFFT properly installed?\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            else:
                logfile.write("Multiple alignments complete\n" \
                       "MAFFT alignments are in files " \
                       "'*.aligned' in "+RESDIR+"\n")

    except fluff.PathError, e:
        logfile.write(e.message+"\n")
    logfile.line()
else:
    logfile.write("Not aligning clusters\n")
    logfile.line()

## ---------------------------------------------------------------------------- ##
#                       |          |                                             #
#                       |          |                                             #
#                       |          |                                             #
#                       |          |                                             #
#                      #            #                                            #
#                      #            #                                            #
#                       ------------                                             #


if options.hmmsearch and options.extracthits and options.blastclust and options.alignment:
    total_residues = fluff.count_nucleotides_from_hmmsearch(args)
    infostring = ("                        PIPELINE STATISTICS\n"+
                  "Total number of characters searched (nucleotides/residues): "+
                  str(total_residues)+"\n")
    logfile.write(infostring)
    logfile.line()

# Cleanup and exit
t = time.asctime(time.localtime())
logfile.write("Cleaning up, pipeline finished at: "+t+"\n")
logfile.line()
fluff.cleanup(TMPDIR)
exit(0)
