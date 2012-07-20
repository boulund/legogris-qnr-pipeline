#!/usr/bin/env python
# Program : qnrpipeline.py
# Author  : Fredrik Boulund 
# Creation date: 2010-07-26
# Release date: 2012-0?-??
# Dependencies:
#   blastclust
#   fluff.py(c)
#   HMMER
#   MAFFT
#   cdbfasta / cdbyank
#
#   This file is the main component of a pipeline that searches databases using 
# HMMER to find sequences matching a hidden Markov model. These sequences are 
# then clustered and aligned against each other to simplify identification of new 
# gene variants and maybe even new genes entirely.
#
# Copyright (C) 2012 Fredrik Boulund
#
# HERE WILL BE LICENSE INFORMATION 

from sys import argv, exit
from os import path, makedirs, system
from datetime import date
from math import ceil, floor
import time
import fluff
import pickle
import shlex, subprocess
from optparse import OptionParser, OptionGroup

ver = "QNR-search pipeline, version 0.8067 BETA" # 2012-07-16
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
try:
    t = time.asctime(time.localtime())
    logfile = open('qnrsearch.log','a')
    if logfile.tell() == 0:
        print "Logfile 'qnrsearch.log' created on",t 
        logfile.write("Logfile 'qnrsearch.log' created on "+t+"\n")
    else:
        print "Logging to 'qnrsearch.log' started on",t
        logfile.write("Logging started on "+t+"\n")
except IOError:
    print "NOTE: cannot create logfile:", arg
    print "Messages will be printed to STDOUT exclusively"

print "\n     -----=====[ ",parser.version," ]=====-----\n\nAssigned work:"
logfile.write("\n     -----=====[  "+parser.version+"  ]=====-----\n\nAssigned work:\n")


## INITIALIZATION, PATH CREATION ETC - - - - - - - - - - - - - - - - - - - - - -
# This is a very big messy pile of if, elif, else statements... As long as it 
# works as intended, I don't feel like refactoring it :).
# Only run HMMsearch. -H
if options.hmmsearch and not(options.extracthits or options.blastclust or options.alignment):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for database(s) to search"
        logfile.write("ERROR! Needs filename(s) for database(s) to search\n")
        parser.print_help()
        exit(1)
    print "Only perform hmmsearch\n"
    logfile.write("Only perform hmmsearch\n")
    options.extracthits = False
    options.blastclust = False
    options.alignment = False
    if not path.isdir(path.abspath(options.hmmsearch_outdir)):
        makedirs(options.hmmsearch_outdir)
        print options.hmmsearch_outdir,"created..."
        logfile.write(options.hmmsearch_outdir+" created...\n")
    else:
        print options.hmmsearch_outdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.hmmsearch_outdir+" already exists, possibly overwriting contents, continuing...\n")

# BOTH run HMMsearch and extract hits, don't cluster or align. -HE
elif (options.hmmsearch and options.extracthits) and not(options.blastclust or options.alignment):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for database(s) to search"
        logfile.write("ERROR! Needs filename(s) for database(s) to search\n")
        parser.print_help()
        exit(1)
    print "Perform hmmsearch and extract hits\n"
    logfile.write("Perform hmmsearch and extract hits\n")
    options.extracthits = True
    options.blastclust = False
    options.alignment = False
    if not path.isdir(path.abspath(options.hmmsearch_outdir)):
        makedirs(options.hmmsearch_outdir)
        print options.hmmsearch_outdir,"created..."
        logfile.write(options.hmmsearch_outdir+" created...\n")
    else:
        print options.hmmsearch_outdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.hmmsearch_outdir+" already exists, possibly overwriting contents, continuing...\n")
    
# Extract hits from hmmsearch AND run clustering, do not align clusters. -EL
elif (options.extracthits and options.blastclust) and not(options.hmmsearch or options.alignment):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for hmmsearch output file(s) to search"
        logfile.write("ERROR! Needs filename(s) for hmmsearch output file(s) to search\n")
        parser.print_help()
        exit(1)
    print "Extract hits from hmmsearch output AND cluster the sequences AND align clusters\n"
    logfile.write("Extract hits from hmmsearch output AND cluster the sequences AND align clusters\n")
    options.hmmsearch = False
    options.alignment = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")
    
# HMMsearch, extract hits and cluster, don't align. -HEL
elif (options.hmmsearch and options.extracthits and options.blastclust) and not(options.alignment):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for database(s) to search"
        logfile.write("ERROR! Needs filename(s) for database(s) to search\n")
        parser.print_help()
        exit(1)
    print "Perform hmmsearch, extract hits and cluster\n"
    logfile.write("Perform hmmsearch, extract hits and cluster\n")
    options.extracthits = True
    options.blastclust = True
    options.alignment = False
    if not path.isdir(path.abspath(options.hmmsearch_outdir)):
        makedirs(options.hmmsearch_outdir)
        print options.hmmsearch_outdir,"created..."
        logfile.write(options.hmmsearch_outdir+" created...\n")
    else:
        print options.hmmsearch_outdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.hmmsearch_outdir+" already exists, possibly overwriting contents, continuing...\n")

# Extract hits from hmmsearch AND run clustering AND align clusters. -ELA
elif (options.extracthits and options.blastclust and options.alignment) and not(options.hmmsearch):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for hmmsearch output file(s) to search"
        logfile.write("ERROR! Needs filename(s) for hmmsearch output file(s) to search\n")
        parser.print_help()
        exit(1)
    print "Extract hits from hmmsearch output AND cluster the sequences\n"
    logfile.write("Extract hits from hmmsearch output AND cluster the sequences\n")
    options.hmmsearch = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")

# Only extract hits from hmmsearch output. -E
elif options.extracthits and not(options.hmmsearch or options.blastclust or options.alignment):
    # Check for filename arguments
    if len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for hmmsearch output file(s) to search"
        logfile.write("ERROR! Needs filename(s) for hmmsearch output file(s) to search\n")
        parser.print_help()
        exit(1)
    print "Only extract hits from hmmsearch output\n"
    logfile.write("Only extract hits from hmmsearch output\n")
    options.hmmsearch = False
    options.blastclust = False
    options.alignment = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")
    
# Only run clustering on pre-existing files. -L
elif options.blastclust and not(options.hmmsearch or options.extracthits or options.alignment):
    print "Only performing clustering\n"
    logfile.write("Only performing clustering\n")
    options.hmmsearch = False
    options.extracthits = False
    options.alignment = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")

# Only run clustering on pre-existing files AND align clusters. -LA
elif (options.blastclust and options.alignment) and not(options.hmmsearch or options.extracthits):
    print "Only performing clustering and cluster alignment\n"
    logfile.write("Only performing clustering and cluster alignment\n")
    options.hmmsearch = False
    options.extracthits = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")

# Only run alignment of pre-existing clusters. -A
elif options.alignment and not(options.hmmsearch or options.extracthits or options.blastclust):
    print "Only performing cluster alignment\n"
    logfile.write("Only performing cluster alignment\n")
    options.hmmsearch = False
    options.extracthits = False
    options.blastclust = False
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")

# If no "only" options were set, run entire pipeline! -HELA
elif (not(options.hmmsearch) and not(options.blastclust) and not(options.extracthits) and not(options.alignment)) or (options.hmmsearch and options.blastclust and options.extracthits and options.alignment):
    # Check for filename arguments
    if args == ['-']:
        print "ERROR! Needs filename(s) for database(s) to search"
        logfile.write("ERROR! Needs filename(s) for database(s) to seach\n")
        parser.print_help()
        exit(1)
    elif len(args) != 0:
        NFILES = len(args)
    else:
        print "ERROR! Needs filename(s) for database(s) to search"
        logfile.write("ERROR! Needs filename(s) for database(s) to search\n")
        parser.print_help()
        exit(1)
    print "Run hmmsearch, parse its output, cluster the results and align clusters (entire pipeline)\n"
    logfile.write("Run hmmsearch, parse its output, cluster the results and align clusters (entire pipeline)\n\n")
    options.hmmsearch = True
    options.blastclust = True
    options.extracthits = True
    options.alignment = True
    if not path.isdir(path.abspath(options.resdir)):
        makedirs(options.resdir)
        print options.resdir,"created..."
        logfile.write(options.resdir+" created...\n")
    else:
        print options.resdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.resdir+" already exists, possibly overwriting contents, continuing...\n")
    if not path.isdir(path.abspath(options.hmmsearch_outdir)):
        makedirs(options.hmmsearch_outdir)
        print options.hmmsearch_outdir,"created..."
        logfile.write(options.hmmsearch_outdir+" created...\n")
    else:
        print options.hmmsearch_outdir,"already exists, possibly overwriting contents, continuing..."
        logfile.write(options.hmmsearch_outdir+" already exists, possibly overwriting contents, continuing...\n")


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
    print TMPDIR,"created..."
    logfile.write(TMPDIR+" created...\n")
else:
    print TMPDIR,"already exists, possibly overwriting contents, continuing..."
    logfile.write(TMPDIR+" already exists, possibly overwriting contents, continuing...\n")

# Make sure that the path to the reference sequences is a complete path
if options.addrefseq:
    if options.addrefseq.startswith("~"):
        QNR_REFERENCE_SEQUENCES_PATH = path.expanduser(options.addrefseq)
    else:
        QNR_REFERENCE_SEQUENCES_PATH = path.abspath(options.addrefseq)


RETR_SEQ_FILEPATH = path.abspath(''.join([TMPDIR,"retrieved_sequences.fasta"]))
RESDIR = path.relpath(options.resdir)

# A minimum score of 0 means that we will consider and classify all hits
# in the hmmsearch output files
MIN_SCORE = options.minscore




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
print settings_string, # string ends with newline so don't print one!

if len(args) > 1:
    printstring = "\nThe following database files were entered at command line:\n"+'\n'.join(args)
else:
    printstring = "\nThe following database file was entered at command line:\n"+'\n'.join(args)

print printstring 
logfile.write(printstring+"\n")
logfile.flush()



logfileseparator = "----------------------------------------------------------------------"
print logfileseparator
logfile.write(logfileseparator+"\n")








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

if options.hmmsearch:
    ##---------------------------------------------------------------------------##
    ## Run hmmsearch on the entered files, assuming they are databases in fasta  ##
    ## format.                                                                   ##
    ##---------------------------------------------------------------------------##
    
    # Path to MODEL:
    model = options.model

    # Path to output directory
    hmmsearch_outdir = options.hmmsearch_outdir

    # CPU flag; "--cpu 4" means four CPUs, empty one core
    cpuflag = ''.join(["--cpu ",str(options.numcpu)])

    # Text wrap long lines; "--notextw" means true, empty false
    # THIS MUST BE ENABLED TO ENSURE CORRECT BEHAVIOR!
    textwflag = "--notextw"

    # Heuristics on/off; --max means no heuristics (max sensitivity), empty full heuristics
    # There is little reason not to use heuristics, HMMer has a higher propensity
    # for crashing if not used and it only increases the number of really low
    # scoring hits anyway...
    if options.noheuristics:
        heurflag = "--max"
    else:
        heurflag = ""

    # Retrieve the path to the model from user set variables above
    model = path.abspath(model) #Using the absolute path (better?)

    # Retrieve input database paths
    # Note the change in usage of variable 'args'! It will soon contain hmmsearch outputfile paths
    databases = args
    args = [] 

    # Retrieve current date, used in output filename to unique:ify the output filenames
    d = date.today()
    t = time.asctime(time.localtime())

    ## HMMsearch with the settings defined above
    # Print a log of the hmmsearch run settings
    logfile.write("Running hmmsearch at:"+t+"\n")
    print "Running hmmsearch at:",t
    #logfile.write("Model used                       : "+path.basename(model)+"\n")
    #logfile.write("Output directory                 : "+hmmsearch_outdir+"\n")
    #logfile.write("CPU-flag (no flag means one cpu) : "+cpuflag+"\n")
    #logfile.write("Text Wrap (empty means textwrap) : "+textwflag+"\n")
    #logfile.write("Sensitivity (empty means default): "+heurflag+"\n")
    #logfile.write("\nThe following input file(s) were entered at command line:\n"+'\n'.join(databases)+"\n")
    #print "Model used                       :", path.basename(model)
    #print "Output directory                 :", hmmsearch_outdir
    #print "CPU-flag (no flag means one cpu) :", cpuflag
    #print "Text Wrap (empty means textwrap) :", textwflag
    #print "Sensitivity (empty means default):", heurflag
    #print "\nThe following input file(s) were entered at command line:\n", '\n'.join(databases)

    for database in databases:
        # Retrieve the filename and path of the database
        outfilename = path.basename(database)
        database = path.abspath(database)

        # Put together the entire string to call hmmsearch 
        call_list = ''.join(["hmmsearch ", cpuflag, " ", textwflag, " ", heurflag, " ",
                             "-o ", hmmsearch_outdir, outfilename, 
                             ''.join([".hmmsearched--", d.isoformat(), " "]), 
                             model, " ", database])
        hmmsearch = shlex.split(call_list)
        # Run hmmsearch
        try:
            output = subprocess.Popen(hmmsearch, stdin=subprocess.PIPE,
                                           stderr=subprocess.PIPE).communicate()
            if "Error: Failed to open hmm file" in output[1]:
                print "CATASTROPHIC: Could not open HMM:",model
                print "Make sure 'model.hmm' is available in current directory or"
                print "supply the -m argument with path to your HMM file"
                logfile.write("CATASTROPHIC: Could not open HMM: "+model+"\n")
                logfile.write("Make sure 'model.hmm' is available in current directory or\n")
                logfile.write("supply the -m argument with path to your HMM file\n")
                exit(1)
            if options.noheuristics:
                args.append(hmmsearch[6]) # 5 contains the output file path, used later
            else:
                args.append(hmmsearch[5]) # 5 contains the output file path, used later
            print "Finished hmmsearch on file",database
            logfile.write("Finished hmmsearch on file "+database+"\n")
        except OSError:
            print "Could not open:", database
            logfile.write("Could not open: "+database+"\n")
        logfile.flush()

    # Output some more details for the log
    t = time.asctime(time.localtime())
    print "Finished hmmsearching the databases at:", t
    print logfileseparator
    logfile.write("Finished hmmsearching the databases at: "+t+"\n")
    logfile.write(logfileseparator+"\n")
    logfile.flush()
else:
    print "Not running hmmsearch"
    logfile.write("Not running hmmsearch\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")


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
    
    t = time.asctime(time.localtime())
    print "Starting to extract and classify hits from hmmsearch output at: "+t
    logfile.write("Starting to extract and classify hits from hmmsearch output at: "+t+"\n")

    # Check that all paths given are valid and that files exists in those locations
    hmmsearch_result_files = []
    for filepath in args: 
        # Note that 'args' might have changed contents from the input arguments,
        # depending of if the entire two first parts of the pipeline were run together.
        if path.isfile(path.abspath(filepath)):
            hmmsearch_result_files.append(path.abspath(filepath))
        else:
            print " ERROR: incorrect path for hmmsearch output file:",filepath
            logfile.write(" ERROR: incorrect path for hmmsearch output file: "+filepath+"\n")

    # Open file for writing the retrieved sequences to semi-temporary file
    try:
        retrseqfile = open(RETR_SEQ_FILEPATH,'w')
    except OSError:
        print " ERROR: Could not open file for writing:", RETR_SEQ_FILEPATH
        logfile.write(" ERROR: Could not open file for writing: "+RETR_SEQ_FILEPATH+"\n")
                

    numerrors = 0
    scores_ids = []
    for hmmsearch_result_file in hmmsearch_result_files:
        # Open file for reading
        try:
            try:
                # Parse the hmmsearch output file -- Extract hits above MIN_SCORE (0)
                print "Parsing",hmmsearch_result_file
                logfile.write("Parsing "+hmmsearch_result_file+"\n")
                parsed = fluff.parse_hmmsearch_output(hmmsearch_result_file,MIN_SCORE)
                score_id_tuples, dbpath = parsed # Unpack parsed information
                scores,dscores,ids = zip(*score_id_tuples) # Unzip the scores/IDs
            except ValueError:
                print "Found no sequences with domain score above or equal", MIN_SCORE,
                print " in file:",hmmsearch_result_file
                logfile.write("Found no sequences with domain score above or equal "+str(MIN_SCORE))
                logfile.write(" in file: "+hmmsearch_result_file+"\n")
                numerrors = numerrors + 1 #score_id_tuples = []
                dbpath = ""
                continue # skip to next file to parse
            except fluff.ParseError as e:
                print e.message
                logfile.write(e.message+"\n")
                numerrors = numerrors +1 #score_id_tuples = []
                dbpath = ""
                continue # skip to next file to parse
             
            if options.retrdb:
                # Retrieve hits directly from their source database, if
                # they classify correctly according to the classification
                # function. The database needs an index file for this, 
                # created using cdbfasta (standard settings), if not available
                # things will go bad.
                try: 
                    sequences, errmessages = fluff.retrieve_sequences_from_db(dbpath, ids, dscores, RETR_SEQ_FILEPATH,
                                                     func=classificationfunction,
                                                     longseqcutoff=options.classifyC,
                                                     longseqdef=options.classifyD)
                    # If there were any error messages, print them and continue.
                    if errmessages:
                        for message in errmessages:
                            print message,
                            logfile.write(message)
                        # If there were any sequences returned, write them.
                        for sequence in sequences:
                            retrseqfile.write(''.join([sequence,"\n"]))
                        print "Retrieved "+str(len(sequences))+" full length sequences from database"
                        logfile.write("Retrieved "+str(len(sequences))+" full length sequences from database\n")
                    else:
                        # Write the identified sequences (full length from databse) 
                        # to disk,they have been classified inside the previous 
                        # function and non-qnr like hit sequences were not extracted.
                        for sequence in sequences:
                            retrseqfile.write(''.join([sequence,"\n"]))
                        print "Retrieved "+str(len(sequences))+" full length sequences from database that contains fragments classified as potential hits"
                        logfile.write("Retrieved "+str(len(sequences))+" full length sequences from databases that contains fragments classified as potential hits\n")
                except ValueError:
                    print "No sequences found in the hmmsearch output file", hmmsearch_result_file
                    logfile.write("No sequences found in the hmmsearch output file "+hmmsearch_result_file+"\n")

            elif int(options.extendleft) or int(options.extendright):
                #print "EXTENDING HITS" # DEBUG
                # Extend the hits to the HMMER model and retrieve "a little more"
                # around the edges of the hits, adds the sequence origin to the
                # fasta headers. Only retrieves sequences that classify correctly
                # according to the classification function.
                try:
                    sequences, errmessages = fluff.extend_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                      ids, MIN_SCORE, dbpath,
                                                                      extendleft=options.extendleft,
                                                                      extendright=options.extendright,
                                                                      func=classificationfunction,
                                                                      longseqcutoff=options.classifyC,
                                                                      longseqdef=options.classifyD)
                    # If there were any error messages, print them and continue.
                    if errmessages:
                        for message in errmessages:
                            print message,
                            logfile.write(message)
                        # Write the original aligned domain sequences that were
                        # retrieved instead of the extended sequences. They have
                        # been classified inside the previous function anyway.
                        for sequence in sequences:
                            retrseqfile.write(''.join([sequence,"\n"]))
                        print "Retrieved "+str(len(sequences))+" sequences from hmmsearch output that classified as potential hits"
                        logfile.write("Retrieved "+str(len(sequences))+" sequences from hmmsearch output that classified as potential hits\n")
                    else:
                        # Write the identified sequences (extended domain hits) to disk,
                        # they have been classified inside the previous function and
                        # non-qnr like sequences have been removed.
                        for sequence in sequences:
                            retrseqfile.write(sequence)
                        print "Retrieved "+str(len(sequences))+" extended sequences from database that classified as potential hits"
                        logfile.write("Retrieved "+str(len(sequences))+" extended sequences from database that classified as potential hits\n")
                except ValueError:
                    print "No sequences found in the hmmsearch output file", hmmsearch_result_file
                    logfile.write("No sequences found in the hmmsearch output file "+hmmsearch_result_file+"\n")
                except fluff.PathError, e:
                    print e.message
                    logfile.write(e.message+"\n")

            else:
                # Retrieve interesting sequences from hmmsearch output
                # and write out to RETR_SEQ_FILEPATH
                # This is where the sequences get their origin attached
                # Classify the hmmsearch hits according to classification function
                # and only write sequences to disk if they are classified as 'true'
                try:
                    sequences = fluff.retrieve_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                        ids, MIN_SCORE, dbpath,
                                                                        func=classificationfunction,
                                                                        longseqcutoff=options.classifyC,
                                                                        longseqdef=options.classifyD)
                    
                    # Write the identified sequences (fragments/domains) to disk,
                    # they have been classified inside the previous function and
                    # non-qnr like sequences have been removed.
                    for sequence in sequences:
                        retrseqfile.write(''.join([sequence,"\n"]))
                    print "Retrieved "+str(len(sequences))+" sequences from hmmsearch output that classified as potential hits"
                    logfile.write("Retrieved "+str(len(sequences))+" sequences from hmmsearch output that classified as potential hits\n")
                except ValueError:
                    print "No sequences found in the hmmsearch output file", hmmsearch_result_file
                    logfile.write("No sequences found in the hmmsearch output file "+hmmsearch_result_file+"\n")
                except fluff.PathError, e:
                    print e.message
                    logfile.write(e.message+"\n")

        except IOError:
            print "Could not open file for reading:", hmmsearch_result_file
            print "Continuing with next file..."
            logfile.write("Could not open file for reading: "+hmmsearch_result_file+"\n")
            logfile.write("Continuing with next file...\n")


        scores_ids.append(score_id_tuples)
        logfile.flush()
    
    if numerrors >= NFILES:
        print "CATASTROPHIC: No input files contains any potential hits!?"
        logfile.write("CATASTROPHIC: No input files contains any potential hits!?\n")
        exit(1)

    # Pickle the scores_ids for disconnected usage in other scripts
    # note that scores_ids is a very messy structure!
    data = scores_ids
    pickle_filename = ''.join([TMPDIR,"pickled.hsseq"])
    outpickle = open(pickle_filename,'wb')
    pickle.dump(data,outpickle)
    outpickle.close()
    retrseqfile.close()


    t = time.asctime(time.localtime())
    print "Finished parsing hmmsearch output files and classifying hits at: "+t
    logfile.write("Finished parsing hmmsearch output files and classifying hits at: "+t+"\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")
    logfile.flush()
else:
    print "Not extracting hmmsearch hits"
    logfile.write("Not extracting hmmsearch hits\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")


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
    ##---------------------------------------------------------------------------##
    ##     Convert sequences to a BLAST database and cluster using blastclust    ##
    ##---------------------------------------------------------------------------##
    # If reference sequences are to be added before clustering, add them
    if options.addrefseq:
        if path.isfile(options.addrefseq):
            append_refseq = "cat "+options.addrefseq+" >> "+RETR_SEQ_FILEPATH
            system(append_refseq)
            print "Added reference sequences from "+options.addrefseq+" to set to cluster"
            logfile.write("Added reference sequences from "+options.addrefseq+" to set to cluster\n")
    
    # Shorten the sequences that are too long, since blastclust seems to have
    # trouble aligning very long sequences (e.g. complete sequence genomes).
    try:
        fluff.limit_sequence_length(RETR_SEQ_FILEPATH,64) # limit to 64 columns of sequence
    except fluff.PathError, e:
        print e.message
        print "\nThe clustering part of the pipeline is dependent of files from previous parts in the "+TMPDIR+" directory" 
        logfile.write(e.message+"\n")
        logfile.write("\nThe clustering part of the pipeline is dependent of files from previous parts in the "+TMPDIR+" directory\n")
        exit(1)

    
    # Uniqueify the sequence IDs; needed for blastclust to cluster them since only
    # the first part of the identifier is parsed and thus sequence IDs risk becoming non-
    # unique. It is also a safeguard against redundant data sets.
    # (TODO: might be possible to use python "set" to remove duplicates, but
    # that would require checking of entire sequences to ensure robustness)
    SeqFilename = TMPDIR+"retrieved_sequences.fasta.shortened"
    UnSeqFilename = TMPDIR+"unique_retrieved_sequences.fasta.shortened"
    try:
        unique_sequences = fluff.uniqueify_seqids(SeqFilename,UnSeqFilename)
    except ValueError:
        print "Could not uniqueify the sequence IDs!"
        logfile.write("Could not uniqueify the sequence IDs!\n")
        fluff.cleanup(TMPDIR)
        exit(1)


    # Run formatdb on the file outputted from uniqueify_seqids and 
    # then run blastclust to cluster results (all in one function)
    t = time.asctime(time.localtime())
    print "Creating temporary database and running blastclust at: ",t
    logfile.write("Creating temporary database and running blastclust at: "+t+"\n")
    logfile.flush()
    
    numcores = options.numcpu                   # Number of CPUs to use, 0 means all
    PercentIdentity = options.percent_identity  # Percent identity threshold, range 3-100
    CovThreshold = options.cov_threshold        # Coverage threshold for blastclust, range 0.1-0.99
    try:
        #print "Running blastclust with the following settings:"
        #print " Percent identity:",PercentIdentity,"\n Coverage Threshold:",CovThreshold
        #print " Number of CPUs:",numcores
        #logfile.write("Running blastclust with the following settings:\n" \
        #              " Percent identity: "+str(PercentIdentity)+"\n Coverage Threshold: "+ \
        #              str(CovThreshold)+"\n Number of CPUs: "+str(numcores)+"\n")
        
        blastclust_return_text, blastclust_output = fluff.run_blastclust(UnSeqFilename,PercentIdentity,CovThreshold,numcores)
        
        if "error" in blastclust_return_text:
            print blastclust_output
            logfile.write(blastclust_output+"\n")
            fluff.cleanup(TMPDIR)
            exit(1)
        elif "ERROR" in blastclust_output[1]:
            print blastclust_output[1]
            logfile.write(blastclust_output[1]+"\n")
            fluff.cleanup(TMPDIR)
            exit(1)
        else:
            print blastclust_return_text
            print blastclust_output[0],
            logfile.write(blastclust_return_text+"\n"+blastclust_output[0])
    
    except fluff.PathError, e:
        print e.message, "\n", UnSeqFilename
        logfile.write(e.message+"\n"+UnSeqFilename+"\n")
        fluff.cleanup(TMPDIR)
        exit(1)


    # Parse output from blastclust (and de-uniqueify sequences IDs -- not needed anymore)
    blastclustoutputfile = UnSeqFilename+".clusters"
    try:
        parsedblastclust = fluff.parse_blastclust(blastclustoutputfile)
    except fluff.PathError, e:
        print e.message, "\n", blastclustoutputfile
        logfile.write(e.message+"\n"+blastclustoutputfile+"\n")
        fluff.cleanup(TMPDIR)
        exit(1)
    except ValueError:
        print "ERROR: Found nothing in blastclust output:", blastclustoutputfile
        logfile.write("ERROR: Found nothing in blastclust output: "+blastclustoutputfile+"\n")
        fluff.cleanup(TMPDIR)
        exit(1)
    
    # Deunique:ify sequence IDs parsed from blastclust output,
    # needed only for writing out cluster scores later on
    clusters = fluff.deuniqueify_seqids(parsedblastclust)
        
    
    # Unpickle the scores_ids (needed for being able to run clustering separately)
    # Really unnecessary to do all the time but does not really matter
    try:
        pkfile = open(''.join([TMPDIR,"pickled.hsseq"]),'rb')
        scores_ids = pickle.load(pkfile)
    except IOError:
        print "Could not read the pickled high-scoring sequences (pickled.hsseq)"
        logfile.write("Could not read the pickled high-scoring sequences (pickled.hsseq)\n")
        fluff.cleanup(TMPDIR)
    except pickle.UnpicklingError:
        print "Could not unpickle pickled.hsseq!"
        logfile.write("Could not unpickle pickled.hsseq!\n")
        fluff.cleanup(TMPDIR)
    
    # Output the identified cluster to files,
    # one with clean clusters and one with scores
    clusterfilename = RESDIR+"/identified_clusters"
    withscoresfilename = RESDIR+"/identified_clusters.scores"
    print "The identified clusters are written to:",clusterfilename
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
    print "Finished clustering sequences at: "+t
    logfile.write("Finished clustering sequences at: "+t+"\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")
    logfile.flush()


    # pickle the clusters and scores_ids for others scripts to use as a convenience
    data = (parsedblastclust,scores_ids)
    pickle_filename = ''.join([TMPDIR,"pickled.clusters"])
    outpickle = open(pickle_filename,'wb')
    pickle.dump(data,outpickle)
    outpickle.close()
else:
    print "Not clustering sequences"
    logfile.write("Not clustering sequences\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")



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
        print "Could not read the pickled clusters (pickled.clusters)"
        logfile.write("Could not read the pickled clusters (pickled.clusters)\n")
        fluff.cleanup(TMPDIR)
    except pickle.UnpicklingError:
        print "Could not unpickle pickled.clusters!"
        logfile.write("Could not unpickle pickled.clusters!\n")
        fluff.cleanup(TMPDIR)
    
    try:
        # If QNR_REFERENCE_SEQUENCES_PATH is invalid, do not align against
        # reference sequences and notify user path was invalid.
        if options.alignseqpath:
            # Log that we started aligning
            t = time.asctime(time.localtime())
            print t+"\nPerforming multiple alignment within identified clusters and "\
                  "cluster vs reference genes"
            logfile.write(t+"\nPerforming multiple alignment within identified clusters and "\
                  "cluster vs reference genes\n")
            logfile.flush()

            # Check that the reference sequence path is valid
            if path.isfile(QNR_REFERENCE_SEQUENCES_PATH):
                retcode = fluff.malign_clusters(parsedblastclust, RESDIR, 
                                                QNR_REFERENCE_SEQUENCES_PATH,
                                                TMPDIR)
            else:
                print ("Filename supplied with -R flag was invalid, NOT aligning"
                       " against any reference sequences.")
                logfile.write("Filename supplied with -R flag was invalid, NOT aligning"
                       " against any reference sequences.\n")
                retcode = fluff.malign_clusters(parsedblastclust, RESDIR, "", TMPDIR)
            if retcode == 1: 
                print "No clusters contained more than one member, nothing to align!"
                logfile.write("No clusters contained more than one member, nothing to align!\n")
            elif retcode == 2:
                print "Could not align sequences using 'mafft'. Is MAFFT properly installed?"
                logfile.write("Could not align sequences using 'mafft'. Is MAFFT properly installed?\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            else:
                print ("Multiple alignments complete\n" \
                       "MAFFT alignments are in files " \
                       "'*.aligned' in "+RESDIR)
                logfile.write("Multiple alignments complete\n" \
                       "MAFFT alignments are in files " \
                       "'*.aligned' in "+RESDIR+"\n")
        else:
            # Log that we started aligning
            t = time.asctime(time.localtime())
            print t+"\nPerforming multiple alignment within identified clusters"
            logfile.write(t+"\nPerforming multiple alignment within identified clusters\n")
            logfile.flush()
            
            # Align the clusters 
            retcode = fluff.malign_clusters(parsedblastclust, RESDIR,
                                            refseqpath="",
                                            workdir=TMPDIR)
            if retcode == 1:
                print "No clusters contained more than one member, nothing to align!"
                logfile.write("No clusters contained more than one member, nothing to align!\n")
            elif retcode == 2:
                print "Could not align sequences using 'mafft'. Is MAFFT properly installed?"
                logfile.write("Could not align sequences using 'mafft'. Is MAFFT properly installed?\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            else:
                print ("Multiple alignments complete\n"\
                       "MAFFT alignments are in files "\
                       "'*.aligned' in "+RESDIR)
                logfile.write("Multiple alignments complete\n" \
                       "MAFFT alignments are in files " \
                       "'*.aligned' in "+RESDIR+"\n")
                                             
    except fluff.PathError, e:
        print e.message
        logfile.write(e.message+"\n")

    print logfileseparator
    logfile.write(logfileseparator+"\n")
else:
    print "Not aligning clusters"
    logfile.write("Not aligning clusters\n")
    print logfileseparator
    logfile.write(logfileseparator+"\n")

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
    print infostring,
    print logfileseparator
    logfile.write(infostring)
    logfile.write(logfileseparator+"\n")

# Cleanup and exit
t = time.asctime(time.localtime())
print "Cleaning up, pipeline finished at: "+t
print logfileseparator
logfile.write("Cleaning up, pipeline finished at: "+t+"\n")
logfile.write(logfileseparator+"\n")
fluff.cleanup(TMPDIR)
exit(0)
