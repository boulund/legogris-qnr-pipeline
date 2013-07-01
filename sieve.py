#!/bin/python
from __future__ import print_function

import leveldb
import time

class Sieve(object):
    def __init__(self, params, logfile):
        self.logfile = logfile
        self.init(params)
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

def run_sieve(sieve, paths, logfile, dbenv):
    (indbpath, infilepath, outdbpath, outfilepath) = paths
    indnadb = inprotdb =  outdnadb =  outprotdb = None
    try:
        if hasattr(sieve, 'indbmode'):
            indnadb = db.DB(dbenv)
            inprotdb = db.DB(dbenv)
            if sieve.indbflags:
                indnadb.set_flags(sieve.indbflags)
                inprotdb.set_flags(sieve.indbflags)
            flag = getattr(sieve, 'indbaccess') if hasattr(sieve, 'indbaccess') else db.DB_RDONLY
            indnadb.open(indbpath, dbname='dna', dbtype=sieve.indbmode, flags=flag | db.DB_THREAD )
            inprotdb.open(indbpath, dbname='prot', dbtype=sieve.indbmode, flags=flag | db.DB_THREAD )
        if hasattr(sieve, 'outdbmode'):
            outdnadb = leveldb.LevelDB(outdbpath+'.dna')
            outprotdb = leveldb.LevelDB(outdbpath+'.prot')
        logfile.writeline('Start: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        return sieve.run(indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath)
    finally:
        if not indnadb is None:
            del indandb
        if not inprotdb is None:
            del inprotdb
        if not outdnadb is None:
            del outdnadb
        if not outprotdb is None:
            del outprotdb
        logfile.writeline('Finish: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        logfile.flush()

def _run_sieves(sieves, dbs, files, logfile, dbenv, startindex=0, endindex=-1):
    if endindex == -1:
        endindex = len(sieves)
    for i in xrange(startindex, endindex):
        if isinstance(sieves[i], tuple):
            (s, params) = sieves[i]
            sieve = s.create(params, logfile)
            inpath = run_sieve(sieve, (dbs[i], files[i], dbs[i+1], files[i+1]), logfile, dbenv)
        elif isinstance(sieves[i], list):
            for s in sieves[i]:
                return _run_sieves([s]+sieves[i::], dbs[i-1::], files[i-1::], logfile, dbenv)

def run_sieves(sieves, dbs, files, logfile, dbpath, startindex=0, endindex=-1):
    dbenv = None
    _run_sieves(sieves, dbs, files, logfile, dbenv, startindex, endindex)
