from __future__ import print_function
from os import path, makedirs, system
import shlex
import subprocess
from datetime import date
import time
import json
from bsddb3 import db
from sieve import Sieve
from parser import Parser

def create(params, logfile):
    return HMMSearch(params, logfile)

class HMMSearch(Sieve):

    def init(self, params):
        self.indbaccess = 0
        self.indbflags = None
        self.indbmode = db.DB_HASH
        self.outdbflags = None
        self.outdbmode = db.DB_RECNO
        self.name = 'HMMer search'
        self.param_names = [
            'model_path',
            ('numcpu', 4),
            # Heuristics on/off; --max means no heuristics (max sensitivity), empty full heuristics
            # There is little reason not to use heuristics, HMMer has a higher propensity
            # for crashing if not used and it only increases the number of really low
            ('use_heuristics', True),
            ('classifyK', 0.7778),
            ('classifyC', 109.64),
            ('classifyM', -7.89),
            ('classifyD', 150.64),
            ('minscore', 0)
        ]


    def run(self, indb, infilepath, outdb, outfilepath):

        classificationfunction = lambda L: self.classifyK*L + self.classifyM

        self.hmmsearch(infilepath, outfilepath)

        parser = Parser(self.logfile, classificationfunction)
        result = parser.parse_file(indb, outfilepath, self.minscore, self.classifyC, self.classifyD)
        for sequence in result:
            doc = json.dumps(sequence)
            outdb.append(doc)


    #Runs hmmsearch on FASTA file, store HMMer output in its own format.
    #Needs to be parsed by Parser.
    def hmmsearch(self, infilepath, outfilepath):
        logfile = self.logfile

        heurflag = '--max' if self.use_heuristics else ''
        hmmsearch_outdir = ''

        t = time.asctime(time.localtime())

        logfile.write("Running hmmsearch at:"+t+"\n")
        logfile.debug("Model used                       : "+path.basename(self.model_path)+"\n")
        logfile.debug("Output directory                 : "+hmmsearch_outdir+"\n")
        logfile.debug("Sensitivity (empty means default): "+heurflag+"\n")

        # Put together the entire string to call hmmsearch
        call_list = 'hmmsearch --cpu %d --notextw %s -o %s %s %s' % (self.numcpu, heurflag, outfilepath, self.model_path, infilepath)
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

if __name__ == 'main':
    pass
