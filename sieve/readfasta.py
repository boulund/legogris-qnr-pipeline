#!/usr/bin/python
from __future__ import print_function
import types
import sys
from datetime import datetime
from itertools import takewhile

from parsing.fasta import FASTAParser
from translator import translate_sequence
from sieve import Sieve

def create(params, logfile):
    return FastaReader(params, logfile)

class FastaReader(Sieve):
    def init(self, params):
        self.outdbmode = True
        self.name = 'FASTA translator'
        self.param_names = [('item_limit', 0)]

    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        outfile = open(outfilepath, 'w')
        parser = FASTAParser(self.logfile)
        prots = {}
        dnas = {}
        n = 0
        try:
            for nseq in takewhile(lambda x: not self.item_limit or n <= self.item_limit, parser.parse_fasta(infilepath)):
                n += 1
                id = str(n)
                outdnadb.put(id, nseq['dna'])
                for (frame, dump, out) in translate_sequence(id, nseq['id'], nseq['desc'], nseq['dna']):
                    outprotdb.put(''.join([id, '_', str(frame)]), dump)
                    outfile.write(out)
        finally:
            outfile.close()


if __name__ == "__main__":
    import cProfile
    import berkeley
    #Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
    #TODO: Clean up, get it workin' or somethin'
    def translate_fasta(inpath, outpath):
        outfile = open(outpath, 'w')
        infile = open(inpath,'r')
        outdb = berkeley.open_fragments('n')
        Sieve().run(None, infile, outdb, outfile)
    #inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
    inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'
    cProfile.run("translate_fasta('"+inpath+"', 'test.pfa')")
    #translate_fasta(inpath, 'test.pfa')
