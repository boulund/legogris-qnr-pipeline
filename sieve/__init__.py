#!/bin/python
from __future__ import print_function

import time

class Sieve(object):
    """Abstract base-class for gene-filtering sieve implementations.

    Subclasses must define a member param_names
    """
    def __init__(self, params, logfile, name, param_names):
        self.param_names = param_names
        self.name = name
        self.logfile = logfile
        for param in self.param_names:
            if isinstance(param, str):
                if not param in params:
                    raise Exception('Missing mandatory parameter %s for sieve %s' % (param, self.name))
                param_name = param
            else:
                param_name = param[0]
                default_value = param[1]
            if param_name in params:
                setattr(self, param_name, params[param_name])
            else:
                setattr(self, param_name, default_value)

def run_sieve(sieve, paths, logfile, dbengine):
    (indbpath, infilepath, outdbpath, outfilepath) = paths
    indnadb = inprotdb =  outdnadb =  outprotdb = None
    try:
        if indbpath:
            (indnadb, inprotdb) = dbengine.open(indbpath, truncate=False)
        if outdbpath:
            (outdnadb, outprotdb) = dbengine.open(outdbpath, truncate=True)
        logfile.writeline('Start: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        return sieve.run(indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath)
    finally:
        if not indnadb is None:
            indnadb.close()
        if not inprotdb is None:
            inprotdb.close()
        if not outdnadb is None:
            outdnadb.close()
        if not outprotdb is None:
            outprotdb.close()
        logfile.writeline('Finish: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        logfile.flush()

def _run_sieves(sieves, dbs, files, logfile, dbengine, startindex=0, endindex=-1):
    if endindex == -1:
        endindex = len(sieves)
    for i in xrange(startindex, endindex):
        if isinstance(sieves[i], tuple):
            (s, params) = sieves[i]
            sieve = s.sieve(params, logfile)
            run_sieve(sieve, (dbs[i], files[i], dbs[i+1], files[i+1]), logfile, dbengine)
        elif isinstance(sieves[i], list):
            for j in xrange(len(sieves[i])):
                s = sieves[i][j]
                sfiles = [files[i]] + [f+str(j) for f in files[i+1::]]
                sdbs = [dbs[i]] + [db+str(j) for db in dbs[i+1::]]
                ssieves = [s]+sieves[i+1::]
                _run_sieves(ssieves, sdbs, sfiles, logfile, dbengine)
            return

def run_sieves(sieves, dbs, files, logfile, dbengine, startindex=0, endindex=-1):
    _run_sieves(sieves, dbs, files, logfile, dbengine, startindex, endindex)
