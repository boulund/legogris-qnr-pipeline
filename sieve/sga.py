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
        for
