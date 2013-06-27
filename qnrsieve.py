#!/bin/python
from logger import Logger
import readfasta
import hmmsearch
from sieve import run_sieves

_DEBUG = True

logfile = Logger('qnrsearch.log')
logfile.open()

inpath = 'tutorial/database/ntsmaller_plus_qnr.nfa'
#inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
#inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'

try:
    run_sieves(
        [
            (readfasta, {'item_limit': 0}),
            (hmmsearch, {'model_path': 'tutorial/model.hmm'})
        ],     #sieves
        ['', 'tutorial/fragments.db', 'tutorial/fragments_passed.db'],   #dbs
        [inpath, 'tutorial/test.pfa', 'tutorial/hmmsearch_out'],   #files
        logfile
    )
finally:
    logfile.close()
