from __future__ import print_function
import shlex
import subprocess

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
        ]


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        call_list = 'scripts/sga.sh %s' % infilepath
        call = shlex.split(call_list)
        # Run hmmsearch
        try:
            output = subprocess.Popen(call, stdin=subprocess.PIPE,
                                        stderr=subprocess.PIPE).communicate()
            print(output[1])
            self.logfile.write("Finished SGA alignment on file "+infilepath+"\n")
        except OSError, e:
            self.logfile.write("Could not open one of sga align or "+infilepath+"\n")
            raise e
        self.logfile.flush()
