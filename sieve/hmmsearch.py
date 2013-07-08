from __future__ import print_function
from os import path, makedirs, system
import shlex
import subprocess
from datetime import date
import time
import json
from bsddb3 import db

import translator
from sieve import Sieve
from parsing.hmmer import HMMERParser
from util import sequence_to_fasta

def create(params, logfile):
    return HMMSearch(params, logfile)

class HMMSearch(Sieve):

    def init(self, params):
        self.indbaccess = 0
        self.indbflags = None
        self.indbmode = db.DB_HASH
        self.outdbflags = None
        self.outdbmode = db.DB_HASH
        self.name = 'HMMer search'
        self.param_names = [
            'model_path',
            'hmmsearch_out',
            ('numcpu', 4),
            ('write_only_domain', False), #If output FASTA file should contain entire input sequence or just matching domain.
            # Heuristics on/off; --max means no heuristics (max sensitivity), empty full heuristics
            # There is little reason not to use heuristics, HMMer has a higher propensity
            # for crashing if not used and it only increases the number of really low
            ('use_heuristics', False),
            ('classifyK', 0.7778),
            ('classifyC', 109.64), #the classification cutoff (minimum score) for long sequence. defparam = 75
            ('classifyM', -7.89),
            ('classifyD', 150.64), #the definition for long sequences. defparam = 85
            ('minscore', 0),
            ('minlength', 20), #minimum fragment length allowed.
            ('classificationfunction', lambda L: self.classifyK*L + self.classifyM)
        ]


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):

        self.hmmsearch(infilepath)

        result = HMMERParser(self.logfile).parse_file(inprotdb, self.hmmsearch_out, self.minscore)

        #Put result in db and sequences in outfile
        passed_count = 0
        outfile = open(outfilepath, 'w')
        try:
            for sequence in result:
                if self.classify_sequence(sequence):
                    id = str(sequence.pop('id', None))
                    outprotdb.put(id, json.dumps(sequence))
                    (dnaid, frame) = id.split('_')
                    frame = int(frame)
                    dna = translator.frame_sequence(indnadb.get(dnaid), frame)
                    outdnadb.put(dnaid, dna)
                    if outfilepath.endswith('pfa'):
                        if self.write_only_domain:
                            seq = sequence['protein'][sequence['dstart']:sequence['dfinish']+1]
                        else:
                            seq = sequence['protein']
                    elif outfilepath.endswith('nfa'):
                        if self.write_only_domain:
                            seq = dna[sequence['dstart']*3:(sequence['dfinish']+1)*3]
                        else:
                            seq = dna
                    outfile.write(sequence_to_fasta(id, seq))
                    passed_count += 1
        finally:
            self.logfile.writeline("%d / %d sequences passed the classification function." % (passed_count, len(result)))
            outfile.close()

    #Runs hmmsearch on FASTA file, store HMMer output in its own format.
    #Needs to be parsed by Parser.
    def hmmsearch(self, infilepath):
        logfile = self.logfile

        heurflag = '--max' if self.use_heuristics else ''
        hmmsearch_outdir = ''

        t = time.asctime(time.localtime())

        logfile.write("Running hmmsearch at:"+t+"\n")
        logfile.debug("Model used                       : "+path.basename(self.model_path)+"\n")
        logfile.debug("Output directory                 : "+hmmsearch_outdir+"\n")
        logfile.debug("Sensitivity (empty means default): "+heurflag+"\n")

        # Put together the entire string to call hmmsearch
        call_list = 'hmmsearch --cpu %d --notextw %s -o %s %s %s' % (self.numcpu, heurflag, self.hmmsearch_out, self.model_path, infilepath)
        hmmsearch = shlex.split(call_list)
        # Run hmmsearch
        try:
            output = subprocess.Popen(hmmsearch, stdin=subprocess.PIPE,
                                        stderr=subprocess.PIPE).communicate()
            if "Error: Failed to open hmm file" in output[1] or 'Error: File existence/permissions problem in trying to open HMM file' in output[1]:
                logfile.write("CATASTROPHIC: Could not open HMM: "+self.model_path+"\n")
                exit(1)
            logfile.write("Finished hmmsearch on file "+infilepath+"\n")
        except OSError, e:
            logfile.write("Could not open one of hmmsearch or "+infilepath+"\n")
            raise e
        logfile.flush()

        # Output some more details for the log
        t = time.asctime(time.localtime())
        logfile.write("Finished hmmsearching the databases at: "+t+"\n")
        logfile.line()
        logfile.flush()

    ##---------------------------------------------------------------------------##
    ##                      CLASSIFY SEQUENCE                                    ##
    ##---------------------------------------------------------------------------##
    def classify_sequence(self, sequence):
        """
        Classifies a sequence as Qnr or not.

        Uses the domain_score and a user defined function
        to classify a given sequence as putative Qnr or not.
        Will always return false for fragments sorter that min_length.

        Input::

            sequence_length The sequence to classify

        Returns::

            classification  a boolean determining whether it should be classified
                            recording to the model or not.

        Errors::

            (none)
        """
        sequence_length = len(sequence['protein'])
        longseqcutoff = self.classifyC
        longseqdef = self.classifyD
        # Pretty self-explanatory. Has a range in which the classification
        # function is used.
        return (sequence_length >= longseqdef and sequence['dscore'] >= longseqcutoff) or (sequence_length >= self.minlength and sequence_length < longseqdef and sequence['dscore'] > self.classificationfunction(sequence_length))

    ######################### END classify_qnr

if __name__ == 'main':
    pass
