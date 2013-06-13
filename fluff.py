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


##-----------------------------------------------##
##            UNIQUE:IFY SEQUENCE IDS            ##
##-----------------------------------------------##
def uniqueify_seqids(sequences,filename):
    '''
    Unique:ify sequences identifiers in fasta formatted
    sequences by adding an integer right after the sequence
    identifer at the first space.

    Input::

        sequences   filename with sequences in fasta format to unique:ify
        filename    the output filename for unique sequences

    Returns::

        (None)      The list of strings in fasta format with all sequences
                    inputted sequences, now with added [integer] right after
                    the sequence identifier, is written out to a file for
                    further usage in the OS environment (i.e. in formatdb)

    Errors::

        ValueError  raised if the input filenames are not valid
        PathError   raised if the input filename path is invalid

    '''

    from os import path
    import re

    if not(isinstance(sequences,str)):
        raise ValueError("Input sequences filename to uniqueify_seqids was not a string!")
    if not(isinstance(filename,str)):
        raise ValueError("Output filename to uniqueify_seqids was not a string!")

    # Open file for reading
    try:
        infile = open(path.abspath(sequences),'r')
    except OSError:
        raise PathError("ERROR: Path to file with sequence IDs to unique:ify is invalid!")
    try:
        outfile = open(path.abspath(filename),'w')
    except OSError:
        raise PathError("ERROR: The output filename for unique sequences could not be opened")



    # Compile a regular expression to match the sequence identifier string and
    # put it in two groups, first with sequence id and second with the
    # rest of the sequence identifier line. The unique number will be inserted
    # between the two.
    ss = re.compile(r'^>([\S\.\-\_]+)(.*\n)')  # Removed a \s in (\s.*\n) at the end 201008xx
    number = 1
    for line in infile:
        if line.startswith(">"):
            hit = re.match(ss,line)
            if hit is not None:
                # Append the number and attach the rest of the line
                # and write to outfile
                outfile.write(''.join([">",hit.group(1),"--",str(number),hit.group(2)]))
                number = number + 1
        else:
            # This would be the entire sequence on following lines until next ">"
            outfile.write(line)

    return
############## END uniqueify_seqids




##-----------------------------------------------##
##                 RUN BLASTCLUST                ##
##-----------------------------------------------##
def run_blastclust(filename,PercentIdentity=90,CovThreshold=50,numcores=4):
    '''
    Run formatdb to create a BLAST database, then run
    blastclust on that database to cluster all hits.

    Input::

        filename        filename with sequences with unique ids to cluster.
        PercentIdentity  the percent identity to cluster with, default is 90 %.
        CovThreshold    the minimum length coverage threshold for clustering,
                        default is 50 %
        numcores        optional argument specifying the number of cores to run
                        blastclust on, default is 4 and 0 means all available.

    Returns::

        (None)      Writes output to a file, 'filename.clusters' that contains
                    all identified clusters on each row.

    Errors::

        PathError   raied if there is something wrong with the paths to output
                    or input files.
        ValueError  raised if there is something wrong

    '''

    from os import path
    import shlex, subprocess

    if not(path.exists(path.abspath(filename))):
        raise PathError("ERROR: The path to unique sequence id hit file is incorrect")

    # Create string for calling the subprocess 'formatdb'
    formatdb_call = ''.join(["formatdb -t SignificantHits_UniqueSeqIDs -i ",filename])
    try:
        formatdb = shlex.split(formatdb_call)
        subprocess.call(formatdb)
        return_text = "Created BLAST database with unique sequence IDs"
    except OSError:
        return_text = "The formatdb command could not be run:", formatdb_call

    # Run blastclust in a similar way
    if not isinstance(CovThreshold,str):
        CovThreshold = str(CovThreshold)
    if not isinstance(PercentIdentity,str):
        PercentIdentity = str(PercentIdentity)
    if not isinstance(numcores,str):
        numcores = str(numcores)

    blastclust_call = ''.join(["blastclust -d ",path.abspath(filename),
                              " -S ",PercentIdentity,
                              " -a ",numcores,
                              " -L ",CovThreshold,
                              " -o ",filename,".clusters"])
    try:
        blastclust = shlex.split(blastclust_call)
        blastclust_output = subprocess.Popen(blastclust,\
                                stdout=subprocess.PIPE,\
                                stderr=subprocess.PIPE).communicate()
    except OSError:
        returnstring = ''.join(["ERROR: Blastclust could not be run. Is it properly installed?\n",
                                "The following call was made: $", blastclust_call])
        return ("error",returnstring)


    return (return_text,blastclust_output)
############## END run_blastclust




##-----------------------------------------------##
##          PARSE OUTPUT FROM BLASTCLUST         ##
##-----------------------------------------------##
def parse_blastclust(filename):
    '''
    Parses blastclust output into a nested list structure

    Input::

        filename    filename of blastclust output

    Returns::

        sequenceIDs list of sequence IDs with unique identifiers right after
                    the '>' symbol

    Errors::

        PathError   raised if the file does not exists
        ValueError  rasied if no unique identifiers could be found and removed

    '''

    from os import path

    if not(path.isfile(path.abspath(filename))):
        raise PathError("ERROR: The blastclust output file could not be opened")

    try:
        filepath = path.abspath(filename)
        file = open(filepath,'r')
    except OSError:
        print "ERROR: Could not open", filepath
        exit(1)

    sequenceIDs = []

    for line in file:
        sequenceIDs.append(line.rstrip("\n ").split(" "))


    if sequenceIDs == []:
        raise ValueError

    return sequenceIDs
############## END parse_blastclust





##-----------------------------------------------##
##            DE-UNIQUEIFY SEQUENCE IDS          ##
##-----------------------------------------------##
def deuniqueify_seqids(sequenceIDs,readfile=False):
    '''
    Removes unique identifiers appended to sequence IDs in
    fasta format by 'uniqueify_seqids'.

    Second option is to read a file and de-unique:ify the
    sequence identifiers found in that file;
    the file formats understood are regular fasta and
    clustalw output.

    Input::

        sequenceIDs list of sequence IDs (parsed from blastclust
                    output)
        readfile    boolean determining whether sequenceIDs contains
                    a filename to read and correct rather than
                    a list of sequence IDs parsed from blastclust.

    Returns::

        sequenceIDs list of sequence IDs without unique identifiers
                    right after the '>' symbol
        None        if readfile was true nothing is return, instead
                    changes are written directly to file.

    Errors::

        ValueError  rasied if no unique identifiers could be found
                    and removed
        PathError   raised if there was some error regarding
                    the file to be read.
    '''
    from os import path
    import re

    if readfile:
        try:
            filein = open(path.abspath(sequenceIDs),'r')
            fileout = open(path.abspath(sequenceIDs+".restored"),'w')
        except OSError:
            raise PathError("Could not open multiple alignment file to de-uniqueify!")

        line = filein.readline()

        # FASTA FORMAT
        if line.startswith(">"):
            print "FASTA format detected"
            regex = r'^(>[\S\.-_]+)--\d+(\s.*\n)'
            regex = re.compile(regex)
            hit = re.match(regex,line)
            if hit is not None:
                fileout.write(''.join([hit.group(1), hit.group(2)]))
            else:
                print "Parse error of first line!!"
                return
            while filein:
                line = filein.readline()
                if line == "":
                    break

                hit = re.match(regex,line)
                if hit is not None:
                    fileout.write(''.join([hit.group(1), hit.group(2)]))
                else:
                    fileout.write(line)

        # CLUSTAL FORMAT
        elif line.startswith("CLUSTAL"):
            print "ClustalW format detected"
            regex = r'([\S\-_]+)(--\d+)(\s+.*\n)'
            regex = re.compile(regex)
            fileout.write(line)
            spacesize = 0
            while filein:
                line = filein.readline()
                if line == "":
                    break
                hit = re.match(regex,line)
                if hit is not None:
                    fileout.write(''.join([hit.group(1), hit.group(3)]))
                    spacesize = len(hit.group(2))
                else:
                    fileout.write(line[spacesize:])

        # UNKNOWN FORMAT; DIE
        else:
            print "No known format detected; accepts only CLUSTAL or fasta"

        return

    # NOT READING FROM FILE BUT FROM LIST OF SEQUENCES
    else:
        clean_sequenceIDs = []

        for cluster in sequenceIDs:
            temp_cluster = [] # Temporary storage of cluster members
            for seqID in cluster:
                # Match the number at the end and remove it
                pattern = r'([\S\-_]+)--\d+'
                found = re.match(pattern,seqID)
                if found is not(None):
                    temp_cluster.append(found.group(1))
            clean_sequenceIDs.append(temp_cluster)

        if clean_sequenceIDs == []:
            raise ValueError

        return clean_sequenceIDs
############## END deuniqueify_seqids





##-----------------------------------------------##
## MULTIPLE ALIGNMENT WITHIN AND BETWEEN CLUSTERS##
##-----------------------------------------------##
def malign_clusters(clusters,resdir,refseqpath,workdir):
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


    # USER EDITABLE LIST WITH REFERENCE GENES
    references = ["QnrA","QnrB","QnrC","QnrD","QnrS"] # Can be modified!



    # Read the file with sequences, fixed filename
    # Each element in the list is a sequence
    sequences = []
    tempseqid = ""
    tempseq = ""
    try:
        seqfile = open(''.join([workdir,"unique_retrieved_sequences.fasta.shortened"]),'r')
        for line in seqfile:
            if line.startswith(">"):
                sequences.append(''.join([tempseqid,tempseq]))
                tempseqid = line
                tempseq = ""
            elif not line.startswith(">"):
                tempseq = ''.join([tempseq,line])
        sequences.append(''.join([tempseqid,tempseq]))
    except OSError:
        raise PathError(''.join(["Could not open/find ",workdir,"'unique_retrieved_sequences.fasta'"]))


    ## Not using clustalw any more, the log-file business has been commented out on purpose!
    # Perform multiple alignment between members of clusters and write log
    # of clustalw stuff to 'workdir'/clustalw.log
    #clustalwlogfile = open(path.abspath(''.join([workdir,'clustalw.log'])),'w')

    clusterno = 0
    for cluster in clusters:
        # Perform multiple alignment within cluster if there are
        # more than one members in cluster
        if len(cluster) > 1:
            # Create file and open it for writing with proper name
            clusterno = clusterno + 1
            malignfilepath = path.abspath(''.join([resdir,"/cluster",
                                                 str(clusterno)]))
            malignfile = open(malignfilepath,"w")


            # Find sequences for all sequence IDs in cluster and write to file
            for seqid in cluster:
                for sequence in sequences:
                    if seqid in sequence:
                        malignfile.write(sequence)
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
##              LIMIT SEQUENCE LENGTH            ##
##-----------------------------------------------##
def limit_sequence_length(sequencefile,MAX_LINES=64):
    '''
    Takes a sequence file and shortens all sequences
    it to the MAX_LINES supplied.

    It is divided by
    inserting the sequence ID header ">seqid...."
    after the maximum sequence length position, thus
    creating several smaller sequence segments with
    the same sequence ID and header line.
    Recommended MAX_LINES is 64 (64 lines of 80
    amino acid residues = 5120).
    As the function splits at "\n" characters it
    requires the fasta file to keep sequences on
    several lines - fasta files with sequences on
    a single line will NOT work.

    Input::

        sequencefile    filename string.
        MAX_LINES       maximum number of lines for one sequence, defaults
                        to 64.

    Returns::

        filename (String)   Writes directly to disk to 'sequencefile.shortened'. filename is the filename of the output file.

    Errors::

        PathError   raised if there is something wrong with path

    '''

    from os import path

    if not path.isfile(path.abspath(sequencefile)):
        raise PathError(''.join(["ERROR: Path ",sequencefile," is not a valid file!"]))
    outfilename = sequencefile+'.shortened'

    try:
        file = open(path.abspath(sequencefile),'r')
        outfile = open(path.abspath(outfilename),'w')
    except OSError:
        raise PathError("Could not open sequence and/or output file!")

    # Reads the ENTIRE file into memory - it is assumed
    # that the number of sequences retrieved does not
    # require more memory than what can fit into memory.
    entire_file = file.read()

    # Splits at all sequences but the first (it does not start with a newline)
    # This method introduces an empty line between the first and following
    # sequences if it is split up, assumed non dangerous but still noted.
    sequences = entire_file.split("\n>")
    first = True # To give special treatment to first sequence
    seqout = []
    for sequence in sequences:
        if not sequence == "":
            if first:
                # Already has initial ">"
                lines = sequence.splitlines(True)
                seqout = []
                counter = 0
            else:
                # Add the'\n>' that was removed
                sequence = ''.join(["\n>",sequence])
                lines = sequence.splitlines(True)
                seqout = []
                counter = 0 # Reset, starting with new sequence

            for line in lines:
                # Add lines of sequence until MAX_LINES,
                # then add sequence identifier again and
                # continue until end of sequence
                if counter > MAX_LINES:
                    if first:
                        if not(lines[0].endswith("\n")):
                            seqout.append(''.join([lines[0],"\n"]))
                        else:
                            seqout.append(lines[0])
                    else:
                        # The following sequences
                        # got split at the initial '\n'
                        # an the identifier is thus in
                        # element 1 instead of 0
                        if not(lines[1].endswith("\n")):
                            seqout.append(''.join([lines[1],"\n"]))
                        else:
                            seqout.append(lines[1])
                    if not(line.endswith("\n")):
                        seqout.append(''.join([line,"\n"]))
                    else:
                        seqout.append(line)
                    counter = 0
                else:
                    seqout.append(line)
                    counter = counter + 1
            #seqout.append("") # Probably not needed...
            first = False # Now the first has been processed

            outfile.write(''.join(seqout))

    return outfilename
############## END limit_sequence_length







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
