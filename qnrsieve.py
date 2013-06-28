#!/bin/python
from logger import Logger
import readfasta
import hmmsearch
import blastclust
from sieve import run_sieves

logfile = Logger('qnrsearch.log')
logfile.open()

inpath = 'tutorial/database/ntsmaller_plus_qnr.nfa'
#inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
#inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'

try:
    run_sieves(
        [
            (readfasta, {'item_limit': 0}),
            (hmmsearch, {'model_path': 'tutorial/model.hmm', 'hmmsearch_out': 'tutorial/hmmsearch_out'}),
            (blastclust, {'blastclust_out': 'tutorial/blastclust_out'})
        ],     #sieves
        ['', 'tutorial/fragments.db', 'tutorial/fragments_passed.db', 'tutorial/clusters.db'],   #dbs
        [inpath, 'tutorial/fragments.pfa', 'tutorial/fragments_passed.pfa', 'tutorial/blastclust_in.pfa'],   #files
        logfile,
        0 #DEBUG: Skip first step
    )
finally:
    logfile.close()
