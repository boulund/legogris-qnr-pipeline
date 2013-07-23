from __future__ import print_function
from os import path, makedirs, system
import shlex
import subprocess
from datetime import date
import time
import json
from itertools import ifilter

from util import translator
from sieve import Sieve
from parsing.hmmer import HMMERParser
from util import sequence_to_fasta

class HMMSearch(Sieve):
    """
    Perform HMMER hmmsearch on input file and assigns sequence and domain scores.
    Optionally filters output where only sequences whose domain score pass a given classification function pass.
    """

    def __init__(self, params, logfile):
        """
        Mandatory parameters:
            * model_path (str): Path to .hmm model file to use for hmmsearch.
        Optional parameters:
            * hmmsearch_out(str, 'hmmsearch_out'): Path to temporary hmmsearch output file.
            * numcpu (int, 4): Number of threads to use for hmmsearch
            * write_only_domain (bool, True): If output FASTA file should contain entire input sequence or just matching domain.
            * use_heuristics (bool, True): Heuristics on/off;  there is little reason not to use heuristics, HMMer has a higher propensity for crashing if not used and it only increases the number of really low scoring hits anyway,
            * longseqdef (int, 151): Sequences longer than this will be considered "long" and are thus compared to `longseqcutoff` instead of `classificationfunction`.
            * longseqcutoff (float, 109.64): The classification cutoff (minimum domain score) for long sequences. That is, "long" sequences will only pass if they have a score higher than this.
            * minscore (int, 0): Minimum domain score, regardless of length, to pass a sequence.
            * min_sequence_length (int, 20): Only consider sequences longer than this.
            * max_domain_length (int, 21844): Do not include sequences where the maximum domain is longer than this.
            * max_sequence_length (int, 21844000): Do not include sequences longer than this.
            * classificationfunction (function, lambda L: self.classifyK*L + self.classifyM): Classification function. The domain score of sequences shorter than `longseqdef` are compared to the result of this function, with the sequence length fed as the parameter. Only those with a domain score higher will pass.
            * classifyK (float, 0.7778): Parameter to the default `classificationfunction`.
            * classifyM (float, -7.89): Parameter to the default `classificationfunction`.
        """

        param_names = [
            'model_path',
            ('hmmsearch_out', 'hmmsearch_out'),
            ('numcpu', 4),
            ('write_only_domain', True),
            ('use_heuristics', True),
            ('longseqdef', 151),
            ('longseqcutoff', 109.64),
            ('minscore', 0),
            ('min_sequence_length', 20), #minimum fragment length allowed.
            ('max_domain_length', 21844),
            ('max_sequence_length', 21844000),
            ('classifyK', 0.7778),
            ('classifyM', -7.89),
            ('classificationfunction', lambda L: self.classifyK*L + self.classifyM)
        ]
        Sieve.__init__(self, params, logfile, name='HMMer hmmsearch', param_names=param_names)


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        self.hmmsearch(infilepath)
        parser = HMMERParser(self.logfile)
        seqs = parser.parse_file(inprotdb, self.hmmsearch_out)

        #Put result in db and sequences in outfile
        passed_count = 0
        outfile = open(outfilepath, 'w')
        try:
            for sequence in ifilter(lambda x: x['dscore'] >= self.minscore, seqs):
                if self.classify_sequence(sequence):
                    id = str(sequence.pop('id', None))
                    outprotdb.put(id, json.dumps(sequence))
                    (dnaid, frame) = id.split('_')
                    frame = int(frame)
                    rawdna = indnadb.get(dnaid)
                    outdnadb.put(dnaid, rawdna)
                    if '.pfa' in outfilepath:
                        if self.write_only_domain:
                            seq = sequence['protein'][sequence['dstart']:sequence['dfinish']+1]
                        else:
                            seq = sequence['protein']
                    elif '.nfa' in outfilepath:
                        dna = translator.frame_sequence(rawdna, frame)
                        if self.write_only_domain:
                            seq = dna[(sequence['dstart']-1)*3:(sequence['dfinish'])*3]
                        else:
                            seq = dna
                    outfile.write(sequence_to_fasta(id, seq))
                    passed_count += 1
        finally:
            self.logfile.writeline("%d sequences passed the classification function." % passed_count)
            outfile.close()

    def hmmsearch(self, infilepath):
        """
        Run hmmsearch on FASTA file, store HMMER output in its own format.
        Use parameters passed on sieve initialization.
        Output can be parsed by `parsing.hmmer`.
        """
        logfile = self.logfile

        heurflag = '' if self.use_heuristics else '--max'
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
        output = subprocess.Popen(hmmsearch, stdin=subprocess.PIPE,
                                    stderr=subprocess.PIPE).communicate()
        if "Error: Failed to open hmm file" in output[1] or 'Error: File existence/permissions problem in trying to open HMM file' in output[1]:
            logfile.write("CATASTROPHIC: Could not open HMM: "+self.model_path+"\n")
            exit(1)
        logfile.write("Finished hmmsearch on file "+infilepath+"\n")
        logfile.flush()

        # Output some more details for the log
        t = time.asctime(time.localtime())
        logfile.write("Finished hmmsearching the databases at: "+t+"\n")
        logfile.line()
        logfile.flush()

    def classify_sequence(self, sequence):
        """
        Classify a sequence as interesting or not.

        Uses the domain score and a user defined function
        to classify a given sequence as putative or not.
        Will always return false for fragments outside the range of min_sequence_length and max_sequence_length and for domains longer than max_domain_length.

        Args:
            sequence The sequence to classify

        Returns:
            classification (bool)  Whether input should be classified according to the model or not.

        """
        sequence_length = len(sequence['protein'])
        domain_length = sequence['dfinish']-sequence['dstart']+1
        return sequence_length <= self.max_sequence_length and domain_length <= self.max_domain_length and sequence_length >= self.min_sequence_length and \
            ( \
                (sequence_length >= self.longseqdef and sequence['dscore'] >= self.longseqcutoff) or \
                (sequence_length < self.longseqdef and sequence['dscore'] > self.classificationfunction(sequence_length)) \
            )

sieve = HMMSearch

if __name__ == 'main':
    pass
