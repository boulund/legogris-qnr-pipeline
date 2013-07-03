from os import path, makedirs, system
from datetime import date
import time
import pickle
import json

from util import PathError, ParseError

TMPDIR = "./pipeline_data/"

class Parser:
    def __init__(self, logfile):
        self.logfile = logfile

    #Reinserts sequences with hmm score and dscore in indb
    def parse_file(self, indb, infilepath, minscore):
        logfile = self.logfile
        score_id_tuples = []
        # Open file for reading
        sequences = []
        try:
            # Parse the hmmsearch output file -- Extract hits above minscore (0)
            logfile.write("Parsing "+infilepath+"\n")
            parsed = _parse_hmmsearch_output(infilepath,minscore)
            score_id_tuples, dbpath = parsed # Unpack parsed information
            scores,dscores,ids = zip(*score_id_tuples) # Unzip the scores/IDs
            for (score, dscore, id) in score_id_tuples:
                seq = json.loads(indb[id])
                seq['score'] = score
                seq['dscore'] = dscore
                indb[id] = json.dumps(seq)
                seq['id'] = id
                sequences.append(seq)

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
def _parse_hmmsearch_output(filename,min_score=0):
    '''
    Parses and retrieves sequence score, maximal domain score
    and sequence IDs from an hmmsearch output file.

    Input::


        file  filename string to hmmsearch output
        min_score  a float with minimum best domain score threshold

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
            is found with score >= min_score.
    '''

    from os import path
    import re

    score_id_tuples = []

    file = open(filename,"r")

    try:
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

            # Match to find the sequence score, maximum domain score and sequence ID
            # from the list of results
            pattern = r'^.{13,14}?(\d+\.\d).{18,20}\s(\d+\.\d).{17}\s([\w\d\._|-]*)'
            found = re.match(pattern,line) # Find the first match to the pattern
            # If the DOMAIN score is higher than or equal min_score count it as a
            # preliminary hit and store the sequence ID and scores.
            if found is not None:
                sequence_score = float(found.group(1))
                domain_score = float(found.group(2))
                sequence_id = found.group(3)
                if domain_score >= min_score:
                    score_id_tuples.append((sequence_score, domain_score, sequence_id))
                    continue
            if line.startswith(">>") or line.startswith("//"): # Now we stop reading the file!
                break
    finally:
        file.close()

    if score_id_tuples == []:
        raise ValueError
    else:
        returntuple = (score_id_tuples,dbpath)
        return returntuple
############## END  parse_hmmsearch_output

