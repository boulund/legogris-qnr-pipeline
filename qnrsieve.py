#!/usr/bin/python
from logger import Logger
import readfasta
import hmmsearch
import blastclust
from sieve import run_sieves

logfile = Logger('qnrsearch.log')
logfile.open()

_PROFILE = True
def run():
    #inpath = 'tutorial/database/ntsmaller_plus_qnr.nfa'
    inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
    #inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'

    try:
        run_sieves(
            [
                (readfasta, {'item_limit': 0}),
                (hmmsearch, {'model_path': 'tutorial/model.hmm', 'hmmsearch_out': 'tutorial/hmmsearch_out'}),
                (blastclust, {'blastclust_out': 'tutorial/blastclust_out', 'clusters_out_path': 'tutorial/identified_clusters', 'clusters_with_scores_out_path': 'tutorial/identified_clusters.scores'})
            ],     #sieves
            ['', 'fragments.db', 'fragments_passed.db', 'clusters.db'],   #dbs
            [inpath, 'tutorial/fragments.pfa', 'tutorial/fragments_passed.pfa', 'tutorial/blastclust_in.pfa'],   #files
            logfile,
            'tutorial/db',
            0,1
        )
    finally:
        logfile.close()

if _PROFILE:
    import cProfile
    cProfile.run('run()')
else:
    run()
