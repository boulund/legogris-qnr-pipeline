#!/usr/bin/python
import sys
from util.logger import Logger
from sieve import run_sieves, readfasta, hmmsearch, blastclust, sga
from db import kyoto, level

logfile = Logger('qnrsearch.log')
logfile.open()

dir = sys.argv[1]
_PROFILE = False
def run():
    inpath = [dir+'/database/ntsmallest1_plus_qnr.nfa', dir+'/database/ntsmallest2_plus_qnr.nfa', dir+'/database/ntsmallest3_plus_qnr.nfa']
    #inpath = dir+'/database/ntsmallest1_plus_qnr.nfa'
    #inpath = dir+'/database/ntsmaller_plus_qnr.nfa'
    #inpath = dir+'/database/ntsmall_plus_qnr.nfa'
    #inpath = dir+'/database/ntsubset_plus_7_qnr.nfa'
    #inpath = dir+'/database/nt_plus_7_qnr.nfa'
    #inpath = dir+'/database/qnr_fragmented.nfa'
    #inpath = dir+'/database/india2.fastq.gz'
    #inpath = ['/lagring/boulund/johan_bengtsson/indien-scilife2011/3_120228_AD0J14ACXX_JL30_index20_1.fastq.gz', '/lagring/boulund/johan_bengtsson/indien-scilife2011/3_120228_AD0J14ACXX_JL30_index20_2.fastq.gz']

    try:
        run_sieves(
            [
                (readfasta, {'item_limit': 0}),
                (hmmsearch, {'model_path': dir+'/model.hmm', 'hmmsearch_out': dir+'/hmmsearch_out'}),
                #[
                #    (hmmsearch, {'model_path': dir+'/model.hmm', 'hmmsearch_out': dir+'/hmmsearch_out'}),
                #j    (hmmsearch, {'model_path': dir+'/model.hmm', 'hmmsearch_out': dir+'/hmmsearch_out'})
                #],
                #(blastclust, {'blastclust_out': dir+'/blastclust_out', 'clusters_out_path': dir+'/identified_clusters', 'clusters_with_scores_out_path': dir+'/identified_clusters.scores'})
                (sga, {})
            ],     #sieves
            ['', dir+'/fragments.db', dir+'/fragments_passed.db', dir+'/clusters.db'],   #dbs
            [inpath, dir+'/readfasta.pfa', dir+'/fragments_passed.nfa', dir+'/clusters.pfa'],   #files
            logfile,
            level,
            #kyoto,
            0,3
        )
    finally:
        logfile.close()

if _PROFILE:
    import cProfile
    cProfile.run('run()')
else:
    run()
