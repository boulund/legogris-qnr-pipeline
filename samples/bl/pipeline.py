#!/usr/bin/env python2
import os
import sys
import time
sys.path.insert(1, os.path.abspath(sys.path[0] + '../../..'))
from util.logger import Logger
from sieve import run_sieves, dnareader, hmmsearch, blastclust, sga, multirunner
from db import level


_PROFILE = False

bp = '/lagring/boulund/johan_bengtsson/indien-scilife2011/'
#paths = [
#    ('R1', [bp+'8_111221_AC03V3ACXX_JL27_index3_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL27_index3_1.fastq.gz']),
#    ('R2', [bp+'7_111221_AC03V3ACXX_JL28_index4_1.fastq.gz', bp+'7_111221_AC03V3ACXX_JL28_index4_2.fastq.gz']),
#    ('R3', [bp+'8_111221_AC03V3ACXX_JL29_index6_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL29_index6_2.fastq.gz']),
#    ('R4', [bp+'3_120228_AD0J14ACXX_JL30_index20_1.fastq.gz', bp+'3_120228_AD0J14ACXX_JL30_index20_2.fastq.gz']),
#    ('R5', [bp+'8_111221_AC03V3ACXX_JL31_index8_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL31_index8_2.fastq.gz']),
#    ('R6', [bp+'8_111221_AC03V3ACXX_JL32_index9_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL32_index9_2.fastq.gz'])
#]
paths = [
    #('ntsmall', ['/lagring/edstromr/genomes/ntsmall_plus_qnr.nfa']),
    ('SKU', [bp+'8_111221_AC03V3ACXX_JL34_index11_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL34_index11_2.fastq.gz']),
    ('SKN', [bp+'8_111221_AC03V3ACXX_JL35_index12_1.fastq.gz', bp+'8_111221_AC03V3ACXX_JL35_index12_2.fastq.gz'])
]
#paths = [('ntsmall', ['tutorial/database/ntsmallest1_plus_qnr.nfa', 'tutorial/database/ntsmallest2_plus_qnr.nfa'])]
    #('', [bp+'.gz', bp+'.gz']),
#inpath = [dir+'/database/ntsmallest1_plus_qnr.nfa', dir+'/database/ntsmallest2_plus_qnr.nfa', dir+'/database/ntsmallest3_plus_qnr.nfa']
#inpath = dir+'/database/ntsmallest1_plus_qnr.nfa'
#inpath = dir+'/database/ntsmaller_plus_qnr.nfa'
#inpath = dir+'/database/ntsmall_plus_qnr.nfa'
#inpath = dir+'/database/ntsubset_plus_7_qnr.nfa'
#inpath = dir+'/database/nt_plus_7_qnr.nfa'
#inpath = dir+'/database/qnr_fragmented.nfa'
#inpath = dir+'/database/india2.fastq.gz'
#inpath = ['/lagring/boulund/johan_bengtsson/indien-scilife2011/3_120228_AD0J14ACXX_JL30_index20_1.fastq.gz', '/lagring/boulund/johan_bengtsson/indien-scilife2011/3_120228_AD0J14ACXX_JL30_index20_2.fastq.gz']
#inpath = ['/lagring/boulund/johan_bengtsson/indien-scilife2011/8_111221_AC03V3ACXX_JL27_index3_1.fastq.gz', '/lagring/boulund/johan_bengtsson/indien-scilife2011/8_111221_AC03V3ACXX_JL27_index3_1.fastq.gz']
def run(dir, inpath):
    if not os.path.exists(dir):
        os.makedirs(dir)
    logfile = Logger(dir + '/qnrsearch.log')
    logfile.open()
    try:
        run_sieves(
            [
                (dnareader, {'item_limit': 0}),
                (hmmsearch, {'model_path': sys.path[0]+'/hmm_models/joint.hmm', 'hmmsearch_out': dir+'/hmmsearch_out', 'write_only_domain': False}),
                [
                    (hmmsearch, {'model_path': sys.path[0]+'/hmm_models/classA.hmm', 'hmmsearch_out': dir+'/hmmsearch_out', 'write_only_domain': True}),
                    (hmmsearch, {'model_path': sys.path[0]+'/hmm_models/classC.hmm', 'hmmsearch_out': dir+'/hmmsearch_out', 'write_only_domain': True}),
                    (hmmsearch, {'model_path': sys.path[0]+'/hmm_models/classD.hmm', 'hmmsearch_out': dir+'/hmmsearch_out', 'write_only_domain': True})
                ],
                (sga, {'error_rate': 0.01, 'min_assembly_overlap': 30, 'min_merge_overlap': 25, 'resolve_small': 10, 'parse_output': False  })
            ],
            ['', dir+'/fragments.db', dir+'/fragments_passed.db', dir+'/fragments_passed_second.db',dir+'/clusters.db'],   #dbs
            [inpath, dir+'/dnareader.pfa', dir+'/fragments_passed.nfa', dir+'/fragments_passed_second.nfa',dir+'/clusters.nfa'],   #files
            logfile,
            level,
            0,3
        )
    finally:
        logfile.close()

def runs(paths):
    for (dir, inpath) in paths:
        print('========================')
        print('Running %s at %s.' % (dir, time.asctime(time.localtime())))
        print('========================')
        run(dir, inpath)
        print('========================')
        print('Finished %s at %s.' % (dir, time.asctime(time.localtime())))
        print('========================')

if _PROFILE:
    import cProfile
    cProfile.run('runs(paths)')
else:
    runs(paths)
