#!/bin/python
from __future__ import print_function
from bsddb import db
import readfasta
import time

_DEBUG = True

def run_sieve(sieve, paths):
    (indbpath, infilepath, outdbpath, outfilepath) = paths
    indb = infile = outdb = outfile = None
    try:
        if indbpath:
            indb = db.DB()
            if sieve.indbflags:
                indb.set_flags(sieve.indbflags)
            indb.open(indbpath, sieve.indbmode, db.DB_RDONLY) #RDONLY only a guess!
        if infilepath:
            infile = open(infilepath,'r')
        if outdbpath:
            outdb = db.DB()
            if sieve.outdbflags:
                outdb.set_flags(sieve.outdbflags)
            outdb.open(outdbpath, sieve.outdbmode, db.DB_CREATE)
        if outfilepath:
            outfile = open(outfilepath, 'w')
        if _DEBUG:
            print('Start', sieve.name, time.asctime(time.localtime()))
        return sieve.run(indb, infile, outdb, outfile)
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
    finally:
        if indb:
            indb.close()
        if infile:
            infile.close()
        if outdb:
            outdb.close()
        if outfile:
            outfile.close()
        if _DEBUG:
            print('Finish:', sieve.name, time.asctime(time.localtime()))

def run_sieves(sieves, dbs, files):
    for i in xrange(0, len(sieves)):
        if isinstance(sieves[i], tuple):
            (sieve, params) = sieves[i]
            inpath = run_sieve(sieve.Sieve(params), (dbs[i], files[i], dbs[i+1], files[i+1]))
        elif isinstance(sieves[i], list):
            for s in sieves[i]:
                return run_sieves([s]+sieves[i::], dbs[i-1::], files[i-1::])

if _DEBUG:
    inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
    #inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'
    run_sieves(
        [(readfasta, {})],     #sieves
        ['', 'fragments.db'],   #dbs
        [inpath, 'test.pfa'],   #files
    )
