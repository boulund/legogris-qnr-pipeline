from __future__ import print_function
import os
from util.combinations import combinations


def init(self, params):
    self.param_names = [
        'sieve',
        'params'
    ]

    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        (outdir, outfile) = outfilepath.rsplit('/', 1)
        for params in combinations(self.params):
            soutdir = '/'.join([outdir,  '_'.join(str(x[1]) for x in params)])
            if not os.pathexists(soutdir):
                os.makedirs(soutdir)
            outdnadb.truncate()
            outprotdb.truncate()
            self.sieve(indnadb, inprotdb, infilepath, outdnadb, outprotdb, '/'.join[soutdir, outfile])


