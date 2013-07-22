from __future__ import print_function
import os

from util.combinations import combinations
from sieve import Sieve

class MultiRunner(Sieve):
    def init(self, params):
        self.indbmode = True
        self.outdbmode = True
        self.name = 'MultiRunner'
        self.param_names = [
            'sieve',
            'params'
        ]

    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        (outdir, outfile) = outfilepath.rsplit('/', 1)
        for params in combinations(self.params):
            soutdir = '/'.join([outdir,  '_'.join(str(x) for x in params.values())])
            if not os.path.exists(soutdir):
                os.makedirs(soutdir)

            if outdnadb is not None:
                outdnadb.truncate()
            if outprotdb is not None:
                outprotdb.truncate()
            s = self.sieve.sieve.create(params, self.logfile)
            soutfile = '/'.join([soutdir, outfile])
            s.run(indnadb, inprotdb, infilepath, outdnadb, outprotdb, soutfile)

sieve = MultiRunner
