#!/usr/bin/python
from util.logger import Logger
from sieve import run_sieves, readfasta, hmmsearch, blastclust, sga
from db import kyoto, level

logfile = Logger('qnrsearch.log')
logfile.open()

_PROFILE = False
def run():
    #inpath = 'tutorial/database/ntsmaller_plus_qnr.nfa'
    #inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
    #inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'
    #inpath = 'tutorial/database/nt_plus_7_qnr.nfa'
    inpath = 'tutorial/database/india2.fastq.gz'

    try:
        run_sieves(
            [
                (readfasta, {'item_limit': 0}),
                (hmmsearch, {'model_path': 'tutorial/model.hmm', 'hmmsearch_out': 'tutorial/hmmsearch_out'}),
                #(blastclust, {'blastclust_out': 'tutorial/blastclust_out', 'clusters_out_path': 'tutorial/identified_clusters', 'clusters_with_scores_out_path': 'tutorial/identified_clusters.scores'})
                (sga, {})
            ],     #sieves
            ['', 'tutorial/fragments.nt.db', 'tutorial/fragments_passed.db', 'tutorial/clusters.db'],   #dbs
            [inpath, 'tutorial/readfasta.nt.pfa', 'tutorial/fragments_passed.nfa', 'tutorial/clusters.pfa'],   #files
            logfile,
            level,
            #kyoto,
            2,3
        )
    finally:
        logfile.close()

if _PROFILE:
    import cProfile
    cProfile.run('run()')
else:
    run()
