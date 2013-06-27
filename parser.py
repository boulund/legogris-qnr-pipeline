from os import path, makedirs, system
from datetime import date
import time
import pickle
import json

from fluff import PathError, ParseError
import berkeley

TMPDIR = "./pipeline_data/"

class Parser:
    def __init__(self, logfile, classificationfunction):
        self.classificationfunction = classificationfunction
        self.logfile = logfile

    #Reinserts sequences with hmm score and dscore in indb
    def parse_file(self, indb, infilepath, minscore, classifyC, classifyD):
        logfile = self.logfile
        classificationfunction = self.classificationfunction
        score_id_tuples = []
        # Open file for reading
        sequences = []
        try:
            # Parse the hmmsearch output file -- Extract hits above minscore (0)
            logfile.write("Parsing "+infilepath+"\n")
            parsed = _parse_hmmsearch_output(infilepath,minscore)
            score_id_tuples, dbpath = parsed # Unpack parsed information
            scores,dscores,ids = zip(*score_id_tuples) # Unzip the scores/IDs
            errmessages = []
            for (score, dscore, id) in score_id_tuples:
                seq = json.loads(indb[id])
                classification = _classify_qnr(sequence_length=len(seq['protein']),
                                            domain_score=dscore,
                                            func=self.classificationfunction,
                                            longseqcutoff=classifyC,
                                            longseqdef=classifyD)
                if classification:
                    seq['score'] = score
                    seq['dscore'] = dscore
                    indb[id] = json.dumps(seq)
                    sequences.append(seq)
                else:
                    print "Not classified:", classification

            if errmessages:
                for message in errmessages:
                    logfile.write(message)
            logfile.write("Retrieved "+str(len(sequences))+" full length sequences from database\n")
            # Write the identified sequences (fragments/domains) to disk,
            # they have been classified inside the previous function and
            # non-qnr like sequences have been removed.
            return sequences

        except IOError:
            logfile.write("Could not open file for reading: "+infilepath+"\n")
            logfile.write("Continuing with next file...\n")
        except ValueError:
            logfile.write("Found no sequences with domain score above or equal "+str(minscore))
            logfile.write(" in file: "+infilepath+"\n")
        except PathError, e:
            logfile.write(e.message+"\n")
        except ParseError as e:
            logfile.write(e.message+"\n")
        return ([], score_id_tuples)


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

    # Pretty self-explanatory. Has a range in which the classification
    # function is used, determined by the first if-statement
    return (sequence_length >= longseqdef and domain_score >= longseqcutoff) or (sequence_length >= minlength and sequence_length < longseqdef and domain_score > func(float(sequence_length)))

    if (int(sequence_length) >= int(longseqdef)) and (float(domain_score) >= float(longseqcutoff)):
        return True
    elif int(sequence_length) < minlength: # PREVOUSLY HARDCODED MIN FRAGMENT LENGTH 20
        return False
    else:
        return int(sequence_length) < int(longseqdef) and float(domain_score) > func(float(sequence_length))


######################### END classify_qnr

