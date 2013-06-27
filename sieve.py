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

def run_sieve(sieve, paths, logfile):
    (indbpath, infilepath, outdbpath, outfilepath) = paths
    indb = outdb = None
    try:
        if hasattr(sieve, 'indbmode'):
            indb = db.DB()
            if sieve.indbflags:
                indb.set_flags(sieve.indbflags)
            flag = getattr(sieve, 'indbaccess') if hasattr(sieve, 'indbaccess') else db.DB_RDONLY
            print(flag)
            indb.open(indbpath, sieve.indbmode, flag | db.DB_THREAD )
        if hasattr(sieve, 'outdbmode'):
            outdb = db.DB()
            if sieve.outdbflags:
                outdb.set_flags(sieve.outdbflags)
            flag = getattr(sieve, 'outdbaccess') if hasattr(sieve, 'outdbaccess') else db.DB_CREATE
            outdb.open(outdbpath, sieve.outdbmode, flag | db.DB_THREAD)
        logfile.writeline('Start: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        return sieve.run(indb, infilepath, outdb, outfilepath)
    finally:
        if not indb is None:
            indb.close()
        if not outdb is None:
            outdb.close()
        logfile.writeline('Finish: %s at %s' % (sieve.name, time.asctime(time.localtime())))

def run_sieves(sieves, dbs, files, logfile):
    for i in xrange(0, len(sieves)):
        if isinstance(sieves[i], tuple):
            (s, params) = sieves[i]
            sieve = s.create(params, logfile)
            inpath = run_sieve(sieve, (dbs[i], files[i], dbs[i+1], files[i+1]), logfile)
        elif isinstance(sieves[i], list):
            for s in sieves[i]:
                return run_sieves([s]+sieves[i::], dbs[i-1::], files[i-1::], logfile)
