from __future__ import print_function
import os

from util.combinations import combinations
from sieve import Sieve

class MultiRunner(Sieve):
    """
    Run another sieve multiple times on the same input with different sets of parameters.
    Useful for evaluating results of different combinations of parameters.
    """
    def __init__(self, params, logfile):
        """
        Mandatory parameters:
            * sieve (module): Sieve to run
            * params (dict): Parameter set listing the different values for each parameter on the following form::

                {
                    'param1': [1, 2],
                    'param2': ['A', 'B', 'C']
                }

            The supplied sieve will run once for each possible combination of parameters. The preceding example will thus run 2*3=6 times.
        """
        param_names = [
            'sieve',
            'params'
        ]
        Sieve.__init__(self, params, logfile, name='MultiRunner', param_names=param_names)

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
