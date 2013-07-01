#!/bin/python
from __future__ import print_function
from bsddb3 import db
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
            outdnadb = db.DB(dbenv)
            outprotdb = db.DB(dbenv)
            if sieve.outdbflags:
                outdnadb.set_flags(sieve.outdbflags)
                outprotdb.set_flags(sieve.outdbflags)
            flag = getattr(sieve, 'outdbaccess') if hasattr(sieve, 'outdbaccess') else db.DB_CREATE
            outdnadb.open(outdbpath, dbname='dna', dbtype=sieve.outdbmode, flags=flag | db.DB_THREAD | db.DB_TRUNCATE)
            outprotdb.open(outdbpath, dbname='prot', dbtype=sieve.outdbmode, flags=flag | db.DB_THREAD)
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

def _run_sieves(sieves, dbs, files, logfile, dbenv, startindex=0):
    for i in xrange(startindex, len(sieves)):
        if isinstance(sieves[i], tuple):
            (s, params) = sieves[i]
            sieve = s.create(params, logfile)
            inpath = run_sieve(sieve, (dbs[i], files[i], dbs[i+1], files[i+1]), logfile, dbenv)
        elif isinstance(sieves[i], list):
            for s in sieves[i]:
                return _run_sieves([s]+sieves[i::], dbs[i-1::], files[i-1::], logfile, dbenv)

def run_sieves(sieves, dbs, files, logfile, dbpath, startindex=0):
    dbenv = db.DBEnv()
    dbenv.set_cachesize(4, 0, 2)
    dbenv.open(dbpath, db.DB_CREATE | db.DB_INIT_MPOOL | db.DB_THREAD, 0)
    try:
        _run_sieves(sieves, dbs, files, logfile, dbenv, startindex)
    finally:
        dbenv.close()
