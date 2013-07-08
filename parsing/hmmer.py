from os import path, makedirs, system
from itertools import takewhile
from datetime import date
import time
import pickle
import json

from util import PathError, ParseError

TMPDIR = "./pipeline_data/"

class HMMERParser:
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
            seqs, dbpath = parsed # Unpack parsed information
            for pseq in seqs:
                seq = json.loads(indb.get(pseq['id']))
                seq['score'] = pseq['score']
                seq['dscore'] = pseq['dscore']
                seq['dstart'] = pseq['dstart']
                seq['dfinish'] = pseq['dfinish']
                indb.put(pseq['id'], json.dumps(seq))
                seq['id'] = pseq['id']
                sequences.append(seq)
            return sequences

        except IOError:
            logfile.write("Could not open file for reading: "+infilepath+"\n")
            logfile.write("Continuing with next file...\n")
        except PathError as e:
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
    seqs = {}

    file = open(filename,"r")

    dbpattern = re.compile(r"# target sequence database:\s*([\w\d\./-]*)")
    headerpattern = re.compile(r'^.{13,14}?(\d+\.\d).{18,20}\s(\d+\.\d).{17}\s([\w\d\._|-]*)')
    bodypattern = re.compile(r'\s+\d+\s.\s+(\d+\.\d).+..\s+(\d+)\s+(\d+)\s+..\s+(\d+)\s+(\d+)\s+..')
    try:
        line = file.readline()
        if not line.startswith("# hmmsearch ::"):
            raise ParseError("NOTE: File is not a valid hmmsearch output file: "+filename)

        #First get scores, dscores and ids from header
        for line in takewhile(lambda l: not l.startswith("Domain annotation") and not l.startswith("//"), file):
            # Extract the path to the database where to collect the sequence data
            foundpath = dbpattern.match(line)
            if foundpath is not(None):
                dbpath = path.abspath(foundpath.group(1))

            # Match to find the sequence score, maximum domain score and sequence ID
            # from the list of results
            found = headerpattern.match(line) # Find the first match to the pattern
            # If the DOMAIN score is higher than or equal min_score count it as a
            # preliminary hit and store the sequence ID and scores.
            if found is not None:
                seq = {}
                seq['score'] = float(found.group(1))
                #seq['dscore'] = float(found.group(2))
                dscore = float(found.group(2))
                seq['id'] = found.group(3)
                if dscore >= min_score:
                    seqs[seq['id']] = seq
                    continue

        #Second, get alignment info
        for line in file:
            if line.startswith('>>'):
                seq = seqs[line.split(' ')[1].strip()]
                continue
            hit = bodypattern.match(line)
            if hit is not None and seq is not None:
                dscore = float(hit.group(1))
                if not seq.has_key('dscore') or seq['dscore'] < dscore:
                    seq['dscore'] = dscore
                    seq['dstart'] = int(hit.group(2))
                    seq['dfinish'] = int(hit.group(3))
            if line.startswith('\n'):
                seq = None
    finally:
        file.close()

    if seqs == {}:
        raise ValueError
    else:
        return (seqs.values(),dbpath)
############## END  parse_hmmsearch_output

