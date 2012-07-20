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



##-----------------------------------------------##
##           PARSE HMMSEARCH OUTPUT              ##
##-----------------------------------------------##
def parse_hmmsearch_output(filename,MIN_SCORE=0):
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
##            FIX FASTA FORMATTING               ##
##-----------------------------------------------##
def fixfasta(sequences):
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





##---------------------------------------------------------------------------##
##                      CLASSIFY SEQUENCES AS QNR                            ##
##---------------------------------------------------------------------------##
def classify_qnr(sequence_length, domain_score, func="", longseqcutoff=75, longseqdef=85, minlength=20):
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
def parse_sequence_positions_from_hmmsearch(filepath, seqid_list, min_score=0, func="",
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
                        if classify_qnr(sequencelength, domainscore, func, 
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





##---------------------------------------------------------------------------##
##                 EXTEND SEQUENCES FROM HMMSEARCH                           ##
##---------------------------------------------------------------------------##
def extend_sequences_from_hmmsearch(hmmfilepath, seqid_list, min_score, dbpath,
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
    seqpositions = parse_sequence_positions_from_hmmsearch(hmmfilepath,
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



##-----------------------------------------------##
##        RETRIEVE SEQUENCES FROM HMMSEARCH      ##
##-----------------------------------------------##
def retrieve_sequences_from_hmmsearch(filepath, seqid_list, min_score, dbpath,
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
                                classification = classify_qnr(sequence_length=len(domain),
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
    sequences = fixfasta(sequences)

    return sequences
############## END retrieve_sequences_from_hmmsearch





##---------------------------------------------------------------------------##
##                 RETRIEVE SEQUENCES FROM DATABASE                          ##
##---------------------------------------------------------------------------##
def retrieve_sequences_from_db(dbpath, seqid_list, domain_scores, retr_seq_filepath,
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

        (none)      Writes directly to disk to 'sequencefile.shortened'.

    Errors::

        PathError   raised if there is something wrong with path

    '''

    from os import path

    if not path.isfile(path.abspath(sequencefile)):
        raise PathError(''.join(["ERROR: Path ",sequencefile," is not a valid file!"]))

    try:
        file = open(path.abspath(sequencefile),'r')
        outfile = open(path.abspath(sequencefile+'.shortened'),'w')
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
                
    return 
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

# fluff.ParseError       General error when parsing
class ParseError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)
