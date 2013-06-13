from os import path, makedirs, system
from datetime import date
import time
import pickle

from fluff import PathError, ParseError

TMPDIR = "./pipeline_data/"

class Parser:
    def __init__(self, logfile, classificationfunction):
        self.classificationfunction = classificationfunction
        self.logfile = logfile


    def parse_files(self, files, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright):
        logfile = self.logfile
        t = time.asctime(time.localtime())
        logfile.write("Starting to extract and classify hits from hmmsearch output at: "+t+"\n")

        # Check that all paths given are valid and that files exists in those locations
        hmmsearch_result_files = []
        for filepath in files:
            if path.isfile(path.abspath(filepath)):
                hmmsearch_result_files.append(path.abspath(filepath))
            else:
                logfile.write(" ERROR: incorrect path for hmmsearch output file: "+filepath+"\n")

        # Open file for writing the retrieved sequences to semi-temporary file
        try:
            retrseqfile = open(outfile,'w')
        except OSError:
            logfile.write(" ERROR: Could not open file for writing: "+RETR_SEQ_FILEPATH+"\n")


        numerrors = 0
        scores_ids = []
        for hmmsearch_result_file in hmmsearch_result_files:
            sequences, score_id_tuples = self.parse_file(hmmsearch_result_file, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright)
            if sequences:
                for sequence in sequences:
                    retrseqfile.write(''.join([sequence,"\n"]))
            scores_ids.append(score_id_tuples)
            logfile.flush()

        if numerrors >= len(hmmsearch_result_files):
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
        logfile.write("Finished parsing hmmsearch output files and classifying hits at: "+t+"\n")
        logfile.line()
        logfile.flush()

    def parse_file(self, hmmsearch_result_file, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright):
        logfile = self.logfile
        classificationfunction = self.classificationfunction
        # Open file for reading
        try:
            # Parse the hmmsearch output file -- Extract hits above minscore (0)
            logfile.write("Parsing "+hmmsearch_result_file+"\n")
            parsed = _parse_hmmsearch_output(hmmsearch_result_file,minscore)
            score_id_tuples, dbpath = parsed # Unpack parsed information
            scores,dscores,ids = zip(*score_id_tuples) # Unzip the scores/IDs
            errmessages = []
            if retrdb:
                # Retrieve hits directly from their source database, if
                # they classify correctly according to the classification
                # function. The database needs an index file for this,
                # created using cdbfasta (standard settings), if not available
                # things will go bad.
                logfile.write("Retrieving full length sequences from database")
                sequences, errmessages = _retrieve_sequences_from_db(dbpath, ids, dscores, outfile,
                                                func=classificationfunction,
                                                longseqcutoff=classifyC,
                                                longseqdef=classifyD)
            elif int(extendleft) or int(extendright):
                #print "EXTENDING HITS" # DEBUG
                # Extend the hits to the HMMER model and retrieve "a little more"
                # around the edges of the hits, adds the sequence origin to the
                # fasta headers. Only retrieves sequences that classify correctly
                # according to the classification function.
                logfile.write("Retrieving extended sequences from database that classified as potential hits")
                sequences, errmessages = _extend_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                ids, minscore, dbpath,
                                                                extendleft=extendleft,
                                                                extendright=extendright,
                                                                func=classificationfunction,
                                                                longseqcutoff=classifyC,
                                                                longseqdef=classifyD)
            else:
                # Retrieve interesting sequences from hmmsearch output
                # and write out to RETR_SEQ_FILEPATH
                # This is where the sequences get their origin attached
                # Classify the hmmsearch hits according to classification function
                # and only write sequences to disk if they are classified as 'true'
                logfile.write("Retrieving sequences from hmmsearch output that classified as potential hits")
                sequences = _retrieve_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                    ids, minscore, dbpath,
                                                                    func=classificationfunction,
                                                                    longseqcutoff=classifyC,
                                                                    longseqdef=classifyD)

            if errmessages:
                for message in errmessages:
                    logfile.write(message)
            logfile.write("Retrieved "+str(len(sequences))+" full length sequences from database\n")
            # Write the identified sequences (fragments/domains) to disk,
            # they have been classified inside the previous function and
            # non-qnr like sequences have been removed.
            return (sequences, score_id_tuples)

        except IOError:
            logfile.write("Could not open file for reading: "+hmmsearch_result_file+"\n")
            logfile.write("Continuing with next file...\n")
        except ValueError:
            logfile.write("Found no sequences with domain score above or equal "+str(minscore))
            logfile.write(" in file: "+hmmsearch_result_file+"\n")
        except PathError, e:
            logfile.write(e.message+"\n")
        except ParseError as e:
            logfile.write(e.message+"\n")
        return ([], score_id_tuples)

##-----------------------------------------------##
##            FIX FASTA FORMATTING               ##
##-----------------------------------------------##
def _fixfasta(sequences):
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

    from math import ceil

    outsequences = []
    for sequence in sequences:
        splitstring = sequence.split("\n",1)
        number_of_rows = int(ceil(len(splitstring[1]) / 80.0))
        seq = []
        if number_of_rows > 1:
            for row in xrange(1,number_of_rows):
                seq.append(splitstring[1][:80] + "\n")
                splitstring[1] = splitstring[1][80:]
            seq.append(splitstring[1]) # + "\n")
            seq.insert(0,splitstring[0] + "\n")
            seq = ''.join(seq)
        else:
            seq = sequence
        outsequences.append(seq)

    return outsequences
############## END fixfasta


##-----------------------------------------------##
##        RETRIEVE SEQUENCES FROM HMMSEARCH      ##
##-----------------------------------------------##
def _retrieve_sequences_from_hmmsearch(filepath, seqid_list, min_score, dbpath,
                                    func="", longseqcutoff=75, longseqdef=85, minlength=20):
    '''
    Retrieves aligned sequence parts from the hmmsearch output file
    for domains with scores above min_score.

    This function potentially requires a lot of memory!

    Input::

        filepath        path to hmmsearch output file.
        seqid_list      list with sequence IDs to retrieve.
        min_score       the minimum domain score.
        dbpath          path do the database where the sequences were
                        hmmsearched from.
        func            (lambda) function for use in classify_qnr.
        longseqcutoff   longseqcutoff for use in classify_qnr.
        longseqthresh   longseqthresh for use in classify_qnr.
        minlength       minimum length for use in classify_qnr.

    Returns::

        seqsnutts   a list with sequences parsed from the alignment.

    Errors::

        PathError   raised if the hmmsearch output file can not be found at the
                    specified location.
        ValueError  raised if no sequences are found in the hmmsearch output
                    file.

    '''

    from os import path
    from string import maketrans
    import re

    # TODO: Make function check if --notextw was enabled
    #print ("NOTICE: Can only retrieve sequences from hmmsearch output if "
    #       "hmmsearch option --notextw was enabled!")

    seqid_list = list(seqid_list)
    filepath = path.abspath(filepath)
    if not(path.exists(filepath)):
        raise PathError("ERROR: The path specified from which to parse sequence "
                        "'coordinates' is not valid.")


    # Regular expression to match the score of the aligned domain
    regex_scores = r'\s+== domain\s+\d+\s+score:\s+(\d+\.\d)'
    regex_scores = re.compile(regex_scores)

    # Create a list that will hold the retrieved sequences with their
    # identifiers.
    sequences = []

    # Name of the database and source file, to be appended to the sequence identifier
    # Parsed out of the path to the source file, assumed to be stored accordingly:
    # /[...]/databasename/sourcefile.whatever
    # Extract the database name from the source path
    dbname = path.abspath(dbpath).split("/")[-2]
    # Join the dbname with the source file name
    source_name = path.basename(dbpath).split(".")
    dbname = "_".join([dbname,source_name[0]])


    file = open(filepath,"r")
    while file:
        line = file.readline()
        if line == "":
            break # EOF

        for sequenceID in seqid_list:
            if line.startswith(">>") and (sequenceID in line):
                # Append the database name to the sequence identifier
                sequence_identifier_line = ''.join([">",dbname,"_",line[3:]])
                # Step down five lines;
                # they contain nothing of interest
                for a in xrange(0,5):
                    line = file.readline()

                # We are currently on lines where there might be alignments
                alignments = True
                while alignments:
                    scorematch = re.match(regex_scores,line)
                    if scorematch is not None:
                        current_domain_score = scorematch.group(1)
                        if float(current_domain_score) > float(min_score):
                            # Read two more uninteresting lines with only
                            # the aligned part of the sequence to the HMM
                            file.readline() #Contains nothing of interest
                            file.readline() #Contains nothing of interest
                            line = file.readline() # this is the line we want!
                            regex_domain = r'\s+([\w\.|\_-]+)\s+\d+\s([\w\*-]+)\s+\d+[\s\n]+'
                            domainhit = re.search(regex_domain,line)
                            if domainhit is not None:
                                # Use translate to remove dashes
                                # and change to uppercase
                                tr = maketrans("","")
                                domain = domainhit.group(2).upper().translate(tr,"-*")

                                # PERFORM CLASSIFICATION OF DOMAIN HIT
                                classification = _classify_qnr(sequence_length=len(domain),
                                                            domain_score=current_domain_score,
                                                            func=func,
                                                            longseqcutoff=longseqcutoff,
                                                            longseqdef=longseqdef,
                                                            minlength=minlength)
                                # If the current domain sequence is classified as a potential
                                # true hit, append it to the list of sequences
                                if classification:
                                    sequences.append(''.join([sequence_identifier_line, domain]))

                    line = file.readline()
                    # Break when we reach the end of the alignments
                    # for this sequence ID
                    if line.startswith(">>"):
                        alignments = False #break inner while loop
                    elif line == "":
                        break # EOF, breaks while loop

    #This should never happend
    if not sequences:
        raise ValueError # No sequence was found!

    # Format the sequences so that they conform
    # better to FASTA "standard"
    sequences = _fixfasta(sequences)

    return sequences
############## END retrieve_sequences_from_hmmsearch

##---------------------------------------------------------------------------##
##                 RETRIEVE SEQUENCES FROM DATABASE                          ##
##---------------------------------------------------------------------------##
def _retrieve_sequences_from_db(dbpath, seqid_list, domain_scores, retr_seq_filepath,
                            func="", longseqcutoff=75, longseqdef=85):
    '''
    Retrieves complete source sequences to hits in hmmsearch from their
    source database files.

    It also classifies the hits and discards those
    that do not meet the criteria. It prepends the fasta headers with
    source information.

    Input::

        dbpath          path to database FASTA file.
        seqid_list      a list with sequence IDs to retrieve.
        domain_scores   a list with domain scores for the sequences IDs.
        retr_seq_filepath   path to output file.
        func            (lambda) function for use in classify_qnr, if this is
                        "" then classify_qnr will use a hardcoded function.
        longseqcutoff   longseqcutoff for use in classify_qnr.
        longseqdef      long sequence definition for use in classify_qnr

    Returns::

        sequences       a list with sequences
        errmessages     a list with error messages (if any)

    Errors::

        (none)

    '''

    from sys import argv, exit
    import subprocess, shlex
    from os import path

    # List of strings with error message if any
    errmessages = []

    # Make sure the seqid_list is a list
    seqid_list = list(seqid_list)

    # Name of the database and source file, to be appended to the sequence identifier
    # Parsed out of the path to the source file, assumed to be stored accordingly:
    # /[...]/databasename/sourcefile.whatever
    # Extract the database name from the source path
    dbname = path.abspath(dbpath).split("/")[-2]
    # Join the dbname with the source file name
    source_name = path.basename(dbpath).split(".")
    dbname = "_".join([dbname,source_name[0]])

    # Make sure there is a cdbfasta index file, *.cidx:
    if path.isfile(str(dbpath)+".cidx"):
        usecdbyank = True
    else:
        usecdbyank = False
        errmessages.append("ERROR: No index file "+str(dbpath)+".cidx,\n Skipping\n")

    sequences = []
    for sequence in seqid_list:
        if usecdbyank:
            # Create a system call to cdbyank
            cdbyank_call = ''.join(["cdbyank ", str(dbpath), ".cidx ", "-a",
                                    "'", str(sequence), "'"])
            # Parse the arguments for Popen
            cdbyank_args = shlex.split(cdbyank_call)
            #print cdbyank_call #DEBUG
            #print cdbyank_args #DEBUG
            # Run cdbyank through Popen and let it communicate stdout to
            # stdout variable, it will never contain more than one retrieved
            # source sequence (though this could potentially be a VERY
            # large genome, and might crash the program if too large).
            stdout, stderr = subprocess.Popen(cdbyank_args,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE).communicate()
            # If something went wrong, make possible to inform the user later!
            if stderr:
                errmessages.append("ERROR: cdbyank "+stderr+" trying to retrieve "+sequence+"\n")
                continue

            # Append the database file name and "source"/location to the FASTA
            # header before appending it the the return sequence list
            extended_domain = stdout.split("\n")
            extended_domain[0] = ''.join([">",dbname,"_",extended_domain[0][1:]])
            extended_domain = '\n'.join(extended_domain)
            sequences.append(str(extended_domain))

    return (sequences, errmessages)
######################### END retrieve_sequences_from_db

##-----------------------------------------------##
##           PARSE HMMSEARCH OUTPUT              ##
##-----------------------------------------------##
def _parse_hmmsearch_output(filename,MIN_SCORE=0):
    '''
    Parses and retrieves sequence score, maximal domain score
    and sequence IDs from an hmmsearch output file.

    Input::


        file  filename string to hmmsearch output
        MIN_SCORE  a float with minimum best domain score threshold

    Returns::

        returntuple  nested tuple with the following contents;
            score_id_tuples, dbpath,
            where score_id_tuples contains triplets with
            sequence score, domain score, sequence id,
            and dbpath contains a string with the path
            to the database searched by hmmsearch.

    Errors::

        ParseError  raised if the file does not conform with
            hmmsearch output format (i.e. does not begin
            with '# hmmsearch ::').
        ValueError  raised if no sequence in the hmmsearch output
            is found with score >= MIN_SCORE.
    '''

    from os import path
    import re

    MIN_SCORE = float(MIN_SCORE)

    # Initialize the list of sequence IDs
    sequence_ids = []
    sequence_scores = []
    domain_scores = []


    file = open(filename,"r")
    line = file.readline()
    if not line.startswith("# hmmsearch ::"):
        raise ParseError("NOTE: File is not a valid hmmsearch output file: "+filename)

    while file:
        line = file.readline()
        # Extract the path to the database where to collect the sequence data
        dbpattern = r"# target sequence database:\s*([\w\d\./-]*)"
        foundpath = re.match(dbpattern,line)
        if foundpath is not(None):
            dbpath = path.abspath(foundpath.group(1))
            #print dbpath #Troubleshooting


        # Match to find the sequence score, maximum domain score and sequence ID
        # from the list of results
        pattern = r'^.{13,14}?(\d+\.\d).{18,20}\s(\d+\.\d).{17}\s([\w\d\._|-]*)'
        found = re.match(pattern,line) # Find the first match to the pattern
        # If the DOMAIN score is higher than or equal MIN_SCORE count it as a
        # preliminary hit and store the sequence ID and scores.
        if found is not None and float(found.group(2)) >= float(MIN_SCORE):
            sequence_scores.append(found.group(1))
            domain_scores.append(found.group(2))
            # Attach the database name to the sequence id for easier identification
            sequence_ids.append(found.group(3))
            #print found.group(1),found.group(2),found.group(3) #Troubleshooting
        elif line.startswith(">>") or line.startswith("//"): # Now we stop reading the file!
            break

    file.close()
    score_id_tuples = zip(sequence_scores,domain_scores,sequence_ids)

    if score_id_tuples == []:
        raise ValueError
    else:
        #print "Found",len(sequence_ids),"sequences above score",MIN_SCORE
        returntuple = (score_id_tuples,dbpath)
        return returntuple
############## END  parse_hmmsearch_output

##---------------------------------------------------------------------------##
##                 EXTEND SEQUENCES FROM HMMSEARCH                           ##
##---------------------------------------------------------------------------##
def _extend_sequences_from_hmmsearch(hmmfilepath, seqid_list, min_score, dbpath,
                                    extendleft=0, extendright=0, func="",
                                    longseqcutoff=75, longseqdef=85, minlength=20):
    '''
    Retrieves "a little more" than the matching domain hits from the
    hmmsearch output.

    It also classifies the hits and discards those
    that do not meet the criteria. It prepends the fasta headers with
    source information according to source file name and folder name.

    Input::

        hmmfilepath     path to hmmsearch output file.
        seqid_list      a list with sequence IDs to retrieve.
        min_score       a primitive first classification.
        dbpath          system path to the source FASTA file/database.
        extendleft      an integer with the number of amino acids to retrieve
                        off the lefthand edge of the hmmsearch alignment.
        extendright     an integer with the number of amino acids to retrieve
                        off the righthand edge of the hmmsearch alignment.
        func            (lambda) function for use in classify_qnr, if this is
                        "" then classify_qnr will use a predetermined function.
        longseqcutoff   longseqcutoff for use in classify_qnr.
        longseqdef      long sequence definition for use in classify_qnr
        minlength       minimum fragment length in classify_qnr.

    Returns::

        sequences       a list of sequences retrieved
        message         a string with error message (if any)

    Errors::

        (none)
    '''

    from sys import argv, exit
    import subprocess, shlex
    from os import path

    # List of strings with error message if any
    errmessages = []

    # Make sure the seqid_list is a list
    # This will crash the script or something if the input is not
    # already a list or a tuple.
    seqid_list = list(seqid_list)

    # Set the stepleft and stepright parameters,
    try:
        extendleft = int(extendleft)
    except ValueError:
        errmessages.append("NOTE: Could not interpret extend parameters, used extendleft=0")
        extendleft = 0
    try:
        extendright = int(extendright)
    except ValueError:
        errmessages.append("NOTE: Could not interpret extend parameters, used extendright=0")
        extendright = 0

    # Find the matched domain sequence positions from hmmsearch output
    seqpositions = _parse_sequence_positions_from_hmmsearch(hmmfilepath,
                        seqid_list, min_score, func, longseqcutoff, longseqdef, minlength)

    # Name of the database and source file, to be appended to the sequence identifier
    # Parsed out of the path to the source file, assumed to be stored accordingly:
    # /[...]/databasename/sourcefile.whatever
    # Extract the database name from the source path
    dbname = path.abspath(dbpath).split("/")[-2]
    # Join the dbname with the source file name
    source_name = path.basename(dbpath).split(".")
    dbname = "_".join([dbname,source_name[0]])

    # Check if the database index file exists, *.cidx
    if path.isfile(str(dbpath)+".cidx"):
        usecdbyank = True
    else:
        usecdbyank = False
        errmessages.append("NOTE: cannot find "+str(dbpath)+".cidx,\n falling back to retrieving aligned domain from hmmsearch output\n")

    #print dbpath #troubleshooting
    sequences = []
    for sequence in seqpositions.keys():

        # For each tuple of positions (i.e. domain) that matched
        # against the model, retrieve the requested part from
        # the database. At each index position in the array,
        # there is a list of tuples with sequence positions.
        for domain in seqpositions[sequence]:
            # Extract sequences positions
            leftpos, rightpos = domain #seqpositions[domain][0]

            # Adjust from which positions to retrieve
            # sequences from source database. Correct if
            # outside bounds
            leftpos = int(leftpos)-extendleft
            if leftpos < 1:
                leftpos = 1
            # No need to check for stepping outside sequence length to the right,
            # since cdbyank handles this
            rightpos = int(rightpos)+extendright

            if usecdbyank:
                # Create a system call to cdbyank
                cdbyank_call = ''.join(["cdbyank ", str(dbpath), ".cidx ", "-R -a ",
                                        "'", str(sequence), " ", str(leftpos), " ", str(rightpos), "'"])
                # Parse the arguments for Popen
                cdbyank_args = shlex.split(cdbyank_call)
                #print cdbyank_call #DEBUG
                #print cdbyank_args #DEBUG
                # Run cdbyank through Popen and let it communicate stdout to
                # stdout variable, it will never contain more than one extended
                # retrieved sequence.
                stdout, stderr = subprocess.Popen(cdbyank_args,
                                                  stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE).communicate()
                # If something went wrong, inform the user and retrieve only the
                # aligned part directly from hmmsearch output as fall back procedure!
                if stderr:
                    errmessages.append("ERROR: cdbyank "+stderr+" when trying to retrieve "+sequence+"\n")
                # Append the database name and "source"/location to the FASTA
                # header of the retrieved extended sequnece before appending
                # it the the return sequence list
                extended_domain = stdout.split("\n")
                extended_domain[0] = ''.join([">",dbname,"_",extended_domain[0][1:]])
                extended_domain = '\n'.join(extended_domain)
                sequences.append(str(extended_domain))
            else:
                # Retrieve the sequence in question directly from hmmsearch
                # output instead.
                seq=[] # Empty this so we only have one sequence in here
                try:
                    # Sequence needs to be a list even if it's only 1 seq
                    seq = retrieve_sequences_from_hmmsearch(hmmfilepath, [sequence],
                                    0, dbpath, func, longseqcutoff, longseqdef)
                    sequences.append(seq[0]+"\n")
                except ValueError:
                    errmessages.append("ERROR: Sequence "+sequence+" could not be retrieved.\n")


    return (sequences, errmessages)
######################### END extend_sequences_from_hmmsearch

##---------------------------------------------------------------------------##
##                      CLASSIFY SEQUENCES AS QNR                            ##
##---------------------------------------------------------------------------##
def _classify_qnr(sequence_length, domain_score, func="", longseqcutoff=75, longseqdef=85, minlength=20):
    """
    Classifies a sequence as Qnr or not.

    Uses the domain_score and a user defined function
    to classify a given sequence as putative Qnr or not.
    Contains a hardcoded minimum fragment length of 10
    under which the function will unconditionally return false.

    Input::

        sequence_length an integer with the sequence length.
        domain_score    a float with the domain score for this sequence.
        func            an optional function to determine classification.
        longseqcutoff   the classification cutoff (minimum score) for long qnr
                        sequences.
        longseqdef      the definition for long sequences.
        minlength       minimum fragment length allowed.

    Returns::

        classification  a boolean determining whether it should be classified
                        as Qnr or not.

    Errors::

        (none)
    """

    # Define a hardcoded function if none given
    if func == "":
        k = 0.7778
        m = -7.954
        func = lambda L: k*L + m

    # Pretty self-explanatory. Has a range in which the classification
    # function is used, determined by the first if-statement
    if (int(sequence_length) >= int(longseqdef)) and (float(domain_score) >= float(longseqcutoff)):
        return True
    elif int(sequence_length) < minlength: # PREVOUSLY HARDCODED MIN FRAGMENT LENGTH 20
        return False
    elif int(sequence_length) < int(longseqdef):
        if float(domain_score) > func(float(sequence_length)):
            return True
        else:
            return False
    else:
        return False

######################### END classify_qnr




##-----------------------------------------------##
##   PARSE SEQUENCES POSITIONS FROM HMMSEARCH    ##
##-----------------------------------------------##
def _parse_sequence_positions_from_hmmsearch(filepath, seqid_list, min_score=0, func="",
                                            longseqcutoff=75, longseqdef=85, minlength=20):
    '''
    Retrieves aligned sequence positions from the hmmsearch output file
    for domains that can be classified according to the classification function.

    Calls classify_qnr.

    Input:
        filepath    path to hmmsearch output file.

        seqid_list  list with sequence IDs to retrieve.
        min_score   the minimum domain score, default=0.
        func        (lambda) function used to classify fragments.
        longseqcutoff   parameter for the classification function.
        longseqdef  parameter for the classification function.
        minlength   minimum fragment length in classify_qnr.

    Returns::

        seqsnutts   a dictionary with sequenceIDs as keys to list containing
                    tuples with position information for the domain alignment
                    hit.

    Errors::

        PathError   raised if the hmmsearch output file can
                    not be found at the specified location.
        ValueError  raised if no sequences are found in the
                    hmmsearch output file.
    '''

    from os import path
    import re

    #print ("NOTICE: Can only retrieve sequences from hmmsearch output if "
    #       "hmmsearch option --notextw was enabled!")

    filepath = path.abspath(filepath)
    if not(path.exists(filepath)):
        raise PathError("ERROR: The path specified from which to parse sequence "
                        "'coordinates' is not valid.")


    seqid_list = list(seqid_list)
    file = open(filepath,"r")
    sequences = [] #PREALLOC

    # Regular expression to match the line containing the scores and coordinates of the alignment of the
    # domain to the HMM.
    # example line: (the vvv's designate what's interesting)
    #         vvv                                                   vvv     vvv
        #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
    #   1 !  413.3   4.9  1.6e-125    1e-120       3     217 .]      12     226 ..       9     226 .. 0.99
    #        (   )                                                 (  )    (  )       (  )    (   )
    regex = r'\s+\d+\s.\s+(\d+\.\d).+..\s+(\d+)\s+(\d+)\s+..\s+(\d+)\s+(\d+)\s+..'
    regex = re.compile(regex)

    # Create a dictionary that will contain key -coordinates for sequence alignments
    # from the hmmsearch output. If several coordinates exists (multiple
    # domains match) they will be added to the list for the correct
    # key (sequence id)
    seqsnutts = {}
    while file:
        line = file.readline()
        if line == "":
            break # EOF

        for sequenceID in seqid_list:
            if line.startswith(">>") and (sequenceID in line):

                # THE FOLLOWING WAS CORRECT UNTIL I DISCOVERED A SPECIAL
                # 'FEATURE' IN HMMER3 THAT MAKES THE NEXT STATEMENT INVALID.
                ## Step down three lines;
                ## they contain nothing of interest
                #for a in xrange(1,4):
                #    line = file.readline()

                line = file.readline()
                # We are currently on lines where there are alignments
                alignments = True
                while alignments:
                    hit = re.match(regex,line)
                    if hit is not None:
                        domainscore = hit.group(1)
                        leftpos = hit.group(2)
                        rightpos = hit.group(3)
                        # Compute the sequence length, +1 is important!
                        sequencelength = int(rightpos)-int(leftpos)+1
                        # CLASSIFY SEQUENCES USING CLASSIFICATION FUNCTION
                        if _classify_qnr(sequencelength, domainscore, func,
                                        longseqcutoff, longseqdef, minlength):
                            try:
                                seqsnutts[sequenceID].append((hit.group(2),hit.group(3)))
                            except KeyError:
                                # This happens when this is a new sequenceID
                                seqsnutts[sequenceID] = [(hit.group(2),hit.group(3))]
                    line = file.readline()

                    # Break when we reach the end of the alignments
                    # for this sequence ID (an empty line)
                    if line.startswith("\n"):
                        alignments = False #break inner loop

    #This should never happend
    if seqsnutts == {}:
        raise ValueError # No sequence was found!


    return seqsnutts
############## END parse_sequence_positions_from_hmmsearch
