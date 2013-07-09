from __future__ import print_function
import shlex
import subprocess
import json

from parsing.fasta import FASTAParser
from util import sequence_to_fasta
from sieve import Sieve
from util import PathError

def create(params, logfile):
    return SGAAligner(params, logfile)

class SGAAligner(Sieve):
    def init(self, params):
        self.indbmode = True
        self.outdbmode = True
        self.name = 'SGA Aligner'
        self.param_names = [
            ('result_out_path', ''),
            ('parse_output', True)
        ]


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        call_list = 'scripts/sga.sh %s' % infilepath
        call = shlex.split(call_list)
        # Run SGA pipeline
        try:
            output = subprocess.Popen(call, stdin=subprocess.PIPE,
                                        stderr=subprocess.PIPE).communicate()
            print(output[1])
            self.logfile.write("Finished SGA alignment on file "+infilepath+"\n")
        except OSError, e:
            self.logfile.write("Could not open one of sga align or "+infilepath+"\n")
            raise e
        self.logfile.flush()
        if self.parse_output:
            #There are two sets of output files: One with the final aligned contigs and one with the removed duplicates.
            #Step one: Associate aligned contigs with original fragments
            parser = FASTAParser(self.logfile)
            aligned = [seq for seq in parser.parse_fasta(infilepath+'.result-contigs.fa')]
            outfile = open(outfilepath, 'w')
            try:
                for pseq in parser.parse_fasta(infilepath+'.filter.pass.fa'):
                    seq = json.loads(inprotdb.get(pseq['id']))
                    for aseq in aligned:
                        if pseq['dna'] in aseq['dna']:
                            seq['contig'] = aseq['id']
                            outprotdb.put(pseq['id'], json.dumps(seq))
                            outfile.write(sequence_to_fasta(''.join([pseq['id'], '_', seq['name'], '_', seq['contig']]), seq['protein']))
            finally:
                outfile.close()
