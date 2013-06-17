#================================================#
# fluff.py                                       #
# Author: Fredrik Boulund                        #
# Date: 2012-07-20                               #
# Version 96                                     #
#                                                #
#   Contains the following functions:            #
#    parse_hmmsearch_output                      #
#    parse_sequence_positions_from_hmmsearch     #
#    count_nucleotides_from_hmmsearch            #
#    fixfasta                                    #
#    classify_qnr                                #
#    extend_sequences_from_hmmsearch             #
#    retrieve_sequences_from_hmmsearch           #
#    retrieve_sequences_from_db                  #
#    uniqueify_seqids                            #
#    run_blastclust                              #
#    parse_blastclust                            #
#    deuniqueify_seqids                          #
#    malign_clusters                             #
#    limit_sequence_length                       #
#    retrieve_sequences_from_fasta               #
#    cleanup                                     #
#                                                #
#   Contains the following exceptions:           #
#    PathError                                   #
#    ParseError                                  #
#================================================#
# Copyright (C) 2012 Fredrik Boulund
#
# THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


##-------------------------------------------------------------##
##     COUNT NUCLEOTIDES SEARCHED FROM HMMSEARCH OUTPUT        ##
##-------------------------------------------------------------##
def count_nucleotides_from_hmmsearch(hmmsearchfiles):
    '''
    Takes hmmsearch-files and summarizes the number of nucleotides searched.

    Input::

        hmmsearchfiles  a list of hmmsearchfiles (paths)

    Returns::

        total_residues  an integer with the number of nucleotides searched
    '''

    import re

    regex = re.compile(r'\((\d+) residues\)')

    total_residues = 0

    for hmmsearchfile in hmmsearchfiles:
        residues = 0
        for line in open(hmmsearchfile):
            hit =  re.search(regex,line)
            if hit is not None:
                residues += int(hit.group(1))
        total_residues += residues

    return total_residues
############## END count_nucleotides_from_hmmsearch

def fixfasta(sequence):

    from math import ceil

    splitstring = sequence.split("\n",1)
    number_of_rows = int(ceil(len(splitstring[1]) / 80.0))
    seq = []
    if number_of_rows > 1:
        for row in xrange(1,number_of_rows):
            seq.append(splitstring[1][:80] + "\n")
            splitstring[1] = splitstring[1][80:]
        seq.append(splitstring[1]) # + "\n")
        seq.insert(0,splitstring[0] + "\n")
        return ''.join(seq)
    else:
        return sequence
##-----------------------------------------------##
##            FIX FASTA FORMATTING               ##
##-----------------------------------------------##
def fixfastas(sequences):
    '''
    Takes a list of sequences and tries to correct their
    formatting.

    Designed to be used only in the
    the function retrieve_sequences_from_hmmsearch.

    Input::

        sequences   list of sequences, each in one
                    complete string with \n markers
                    between identifier line and sequence.

    Returns::

        outsequences   list of sequences with hopefully
                    better formatting.

    Errors::

        (none)
    '''


    outsequences = []
    for sequence in sequences:
        seq = fixfasta(sequence)
        outsequences.append(seq)

    return outsequences
############## END fixfasta

def fragment_to_fasta(fragment, sequence='AA'):
    if sequence == 'DNA':
        return fixfasta(''.join(['>', fragment['id'], '\n', fragment['dna'], '\n']))
    if sequence == 'AA': #Amino acids
        return fixfasta(''.join(['>', fragment['id'], '\n', fragment['protein'], '\n']))

##-----------------------------------------------##
## MULTIPLE ALIGNMENT WITHIN AND BETWEEN CLUSTERS##
##-----------------------------------------------##
def malign_clusters(clusters,resdir,refseqpath,seqfilepath):
    '''
    Reads clusters and from file reads complete sequences
    for the sequence IDs specified in the clusters.

    Note that it reads sequences from the filename
    'workdir'/retrieved_sequences.fasta.
    The function assumes that MAFFT is installed and
    accessible through the PATH variable or in the current
    directory. If refseqpath is empty the function does
    not align against any reference sequences. Has hardcoded
    the names of the five plasmid mediated Qnr-genes around line
    1029.


    Input::

        clusters    nested list stucture with clusters.
        resdir      directory to output results.
        refseqpath  path to file with reference qnr sequence in fasta format.
                    Can be an empty string if no reference sequences are to be
                    aligned against.
        workdir     working directory with temporary files
                    (unique_retrieved_sequences.fasta.shortened).

    Returns::

        (nothing)   On success writes multiple alignments to file:
                    resdir/cluster*.aligned
        2 (int)     If there is a MAFFT error.

    Errors::

        PathError   raised if retrieved_sequences.fasta could not be
                    opened/found.

    '''

    from os import path, system
    import re
    import shlex, subprocess


    clusterno = 0
    # USER EDITABLE LIST WITH REFERENCE GENES
    references = ["QnrA","QnrB","QnrC","QnrD","QnrS"] # Can be modified!

    for cluster in clusters.values():
        # Perform multiple alignment within cluster if there are
        # more than one members in cluster
        if len(cluster) > 1:
            # Create file and open it for writing with proper name
            clusterno = clusterno + 1
            malignfilepath = path.abspath(''.join([resdir,"/cluster",
                                                 str(clusterno)]))
            malignfile = open(malignfilepath,"w")

            # Find sequences for all sequence IDs in cluster and write to file
            for fragment in cluster:
                malignfile.write(fragment_to_fasta(fragment))
            malignfile.close()

            # Set parameters etc for clustalw and make alignment
            infilepath = malignfilepath
            outfilepath = ''.join([malignfilepath,".aligned"])
            mafft_call = ''.join(["mafft --quiet ",infilepath," > ",outfilepath])
            errors = system(mafft_call)
            if errors != 0:
                return 2 # Error code for mafft error

            # The following code has been commented out because the
            # usage of mafft is preferred over clustalw. It can be
            # reused if one prefers clustalw.
            # It does produce worse alignments, especially over sequences
            # of very different lengths...
            #clustalw_call = ''.join(["clustalw -INFILE=",infilepath,
            #                        " -OUTFILE=",outfilepath])
            #args = shlex.split(mafft_call)
            #subprocess.call(args)
            #clustalwlog = subprocess.Popen(args,stdout=subprocess.PIPE,
            #                               stderr=subprocess.PIPE).communicate()
            #if not len(clustalwlog[0]) == 0:
            #    clustalwlogfile.write(clustalwlog[0])
    if clusterno == 0:
        #print "\n!! --- No clusters contained more than one member (only singletons)\n"
        return 1 # No need to stay here, lets go home to our box on the highway


    # If refseqpath is empty just skip over
    # this alignment against reference sequences.
    if refseqpath is not "":
        # Read the reference sequences from the path specified
        # Each element in the list is a complete sequence
        refsequences = []
        tempseqid = ""
        tempseq = ""
        try:
            seqfile = open(refseqpath,'r')
            for line in seqfile:
                if line.startswith(">"):
                    refsequences.append(''.join([tempseqid,tempseq]))
                    tempseqid = line
                    tempseq = ""
                elif not line.startswith(">"):
                    tempseq = ''.join([tempseq,line])
            refsequences.append(''.join([tempseqid,tempseq]))
        except OSError:
            raise PathError(''.join(["ERROR: cannot open", refseqpath]))


        # Perform multiple alignment against all reference plasmid mediated qnr-genes
        # for each cluster
        clusterno = 0
        for cluster in clusters:
            if len(cluster) > 1:
                clusterno = clusterno + 1

                for reference in references:
                    malignfilepath = path.abspath(''.join([resdir,"ref",reference,
                                                           ".cluster",str(clusterno)]))
                    malignfile = open(malignfilepath,"w")
                    # Find sequences for all sequence IDs in cluster and write to file
                    for seqid in cluster:
                        for sequence in sequences:
                            if seqid in sequence:
                                malignfile.write(sequence)
                    # Also write the reference gene to the same file
                    for refseq in refsequences:
                        if reference in refseq:
                            malignfile.write(refseq)
                    malignfile.close()

                    # Set parameters etc for clustalw and make alignment
                    infilepath = malignfilepath
                    outfilepath = ''.join([malignfilepath,".aligned"])
                    mafft_call = ''.join(["mafft --quiet ",infilepath," > ",outfilepath])
                    errors = system(mafft_call)
                    if errors != 0:
                        return 2 # Error code for mafft error

                    # The following code has been commented out because the
                    # usage of mafft is preferred over clustalw. It can be
                    # reused if one prefers clustalw. Using clustalw gives
                    # less good alignments...
                    #clustalw_call = ''.join(["clustalw -INFILE=",infilepath,
                    #                        " -OUTFILE=",outfilepath])
                    #args = shlex.split(mafft_call)
                    #clustalwlog = subprocess.Popen(args,stdout=subprocess.PIPE,
                    #                               stderr=subprocess.PIPE).communicate()
                    #if not len(clustalwlog[0]) == 0:
                    #    clustalwlogfile.write(clustalwlog[0])
    return
############## END malign_clusters


##-----------------------------------------------##
##          RETRIEVE SEQUENCES FROM FASTA        ##
##-----------------------------------------------##
def retrieve_sequences_from_fasta(fastapath,seqid_list):
    '''
    Retrieves sequence(s) from a fasta file.


    Python implementation so it is a bit slow, not intended
    for retrieving sequences from large databases.

    Input::

        fastapath   string with path to fasta file
        seqid_list  list with sequence IDs to retrieve

    Returns::

        sequences   list of strings in fasta format with all sequences to
                    retrieve, containing linebreaks after sequence ID headers.

    Errors::

        PathError   raised if the supplied path does not exists or something
                    related to that.
        ValueError  raised if no sequences are found in the database.

    '''

    from os import path

    fastapath = path.abspath(fastapath)
    if not(path.isfile(fastapath)):
        raise PathError("ERROR: Path to fasta file is incorrect.")

    seqid_list = list(seqid_list)
    db = open(fastapath,"r")
    sequences = [] #PREALLOC

    #print "Trying to retrieve",len(seqid_list),"sequences from fasta file:", fastapath
    seq_counter = 0
    line = db.readline()
    while db:
        if line.startswith(">"): #Reached sequence descriptor line
            for id in seqid_list:
                if id in line:
                    seq_counter = seq_counter + 1
                    sequencestring = line
                    seq = True
                    while seq:
                        line = db.readline()
                        if not(line.startswith(">")):
                            sequencestring = ''.join([sequencestring,line])
                        elif line.startswith(">"):
                            seq = False
                    sequences.append(sequencestring)
                else:
                    line = db.readline()

            if seqid_list == []:
                break

        elif line == "": # EOF
            db = False

        else:
            line = db.readline()

    if sequences == []:
        raise ValueError # No sequence was found!

    #print "Retrieved",seq_counter,"sequences from file"
    return sequences
############## END retrieve_sequences_from_fasta




##-----------------------------------------------##
##                   CLEANUP                     ##
##-----------------------------------------------##
def cleanup(tmpdir,exitcode=0):
    '''
    Cleans up (moves) error.log and formatdb.log to
    specified tmp directory

    Input::

        tmpdir      name of the temporary directory to which to move the two
                    log files.
        exitcode    the exitcode to give the shell, defaults to 0.

    Returns::

        (none)

    Errors::

        (none)      Will print error messages to stderr if
                    the files does not exist or could not
                    be moved.

    '''
    import shlex, subprocess
    from sys import exit
    args = shlex.split("mv error.log formatdb.log "+tmpdir)
    out = subprocess.Popen(args, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE).communicate()
    if out[1] != "":
        pass
############## END cleanup


##-----------------------------------------------##
##                CUSTOM EXCEPTIONS              ##
##-----------------------------------------------##

# fluff.PathError       General error with path
class PathError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class ParseError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)
