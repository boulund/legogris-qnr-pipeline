from datetime import date
import time
import pickle

from fluff import PathError, cleanup

TMPDIR = "./pipeline_data/"

class BLASTClusterer:
    def __init__(self, logfile):
        self.logfile = logfile

    # Number of CPUs to use, 0 means all
    # Percent identity threshold, range 3-100
    # Coverage threshold for blastclust, range 0.1-0.99
    def run(self, filepath, numcores, percent_identity, cov_threshold):
        logfile = self.logfile
        # Shorten the sequences that are too long, since blastclust seems to have
        # trouble aligning very long sequences (e.g. complete sequence genomes).
        try:
            SeqFilename = _limit_sequence_length(filepath,64) # limit to 64 columns of sequence
        except PathError, e:
            logfile.write(e.message+"\n")
            logfile.write("\nThe clustering part of the pipeline is dependent of files from previous parts in the "+TMPDIR+" directory\n")
            exit(1)


        # Uniqueify the sequence IDs; needed for blastclust to cluster them since only
        # the first part of the identifier is parsed and thus sequence IDs risk becoming non-
        # unique. It is also a safeguard against redundant data sets.
        # (TODO: might be possible to use python "set" to remove duplicates, but
        # that would require checking of entire sequences to ensure robustness)
        # SeqFilename = TMPDIR+"retrieved_sequences.fasta.shortened"
        # UnSeqFilename = TMPDIR+"unique_retrieved_sequences.fasta.shortened"
        UnSeqFilename = SeqFilename+'.unique'
        try:
            unique_sequences = _uniqueify_seqids(SeqFilename,UnSeqFilename)
        except ValueError:
            logfile.write("Could not uniqueify the sequence IDs!\n")
            cleanup(TMPDIR)
            exit(1)


        # Run formatdb on the file outputted from uniqueify_seqids and
        # then run blastclust to cluster results (all in one function)
        t = time.asctime(time.localtime())
        logfile.write("Creating temporary database and running blastclust at: "+t+"\n")
        logfile.flush()

        try:
            #logfile.write("Running blastclust with the following settings:\n" \
            #              " Percent identity: "+str(percent_identity)+"\n Coverage Threshold: "+ \
            #              str(cov_threshold)+"\n Number of CPUs: "+str(numcores)+"\n")

            blastclust_return_text, blastclust_output = _run_blastclust(UnSeqFilename,percent_identity,cov_threshold,numcores)

            if "error" in blastclust_return_text:
                logfile.write(blastclust_output+"\n")
                cleanup(TMPDIR)
                exit(1)
            elif "ERROR" in blastclust_output[1]:
                logfile.write(blastclust_output[1]+"\n")
                cleanup(TMPDIR)
                exit(1)
            else:
                logfile.write(blastclust_return_text+"\n"+blastclust_output[0])

        except PathError, e:
            logfile.write(e.message+"\n"+UnSeqFilename+"\n")
            cleanup(TMPDIR)
            exit(1)


        # Parse output from blastclust (and de-uniqueify sequences IDs -- not needed anymore)
        blastclustoutputfile = UnSeqFilename+".clusters"
        try:
            parsedblastclust = _parse_blastclust(blastclustoutputfile)
        except PathError, e:
            logfile.write(e.message+"\n"+blastclustoutputfile+"\n")
            cleanup(TMPDIR)
            exit(1)
        except ValueError:
            logfile.write("ERROR: Found nothing in blastclust output: "+blastclustoutputfile+"\n")
            cleanup(TMPDIR)
            exit(1)

        # Deunique:ify sequence IDs parsed from blastclust output,
        # needed only for writing out cluster scores later on
        clusters = _deuniqueify_seqids(parsedblastclust)

        # Unpickle the scores_ids (needed for being able to run clustering separately)
        # Really unnecessary to do all the time but does not really matter
        try:
            pkfile = open(''.join([TMPDIR,"pickled.hsseq"]),'rb')
            scores_ids = pickle.load(pkfile)
        except IOError:
            logfile.write("Could not read the pickled high-scoring sequences (pickled.hsseq)\n")
            cleanup(TMPDIR)
        except pickle.UnpicklingError:
            logfile.write("Could not unpickle pickled.hsseq!\n")
            cleanup(TMPDIR)
        return (clusters, parsedblastclust, scores_ids)

##-----------------------------------------------##
##              LIMIT SEQUENCE LENGTH            ##
##-----------------------------------------------##
def _limit_sequence_length(sequencefile,MAX_LINES=64):
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
##            UNIQUE:IFY SEQUENCE IDS            ##
##-----------------------------------------------##
def _uniqueify_seqids(sequences,filename):
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
def _run_blastclust(filename,PercentIdentity=90,CovThreshold=50,numcores=4):
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
def _parse_blastclust(filename):
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
def _deuniqueify_seqids(sequenceIDs,readfile=False):
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

