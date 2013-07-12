#!/usr/bin/python
from __future__ import print_function
import types
import sys
from datetime import datetime
from itertools import takewhile

from parsing.fasta import FASTAParser
from parsing.fastq import FASTQParser
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
        startid = 0
        open(outfilepath, 'w').close()
        if isinstance(infilepath, list):
            for infile in infilepath:
                print('Opening new file %s, n: %d' % (infile, startid))
                startid = self.run_file(indnadb, inprotdb, infile, outdnadb, outprotdb, outfilepath, startid)
        else:
            self.run_file(indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath)

    def run_file(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath, startid=0):
        outfile = open(outfilepath, 'a')
        if '.gz' in infilepath:
            parser = FASTQParser(self.logfile, gzip=True)
        elif '.fastq' in infilepath or '.fq' in infilepath:
            parser = FASTQParser(self.logfile, gzip=False)
        else:
            parser = FASTAParser(self.logfile)
        prots = {}
        dnas = {}
        n = startid
        try:
            for nseq in takewhile(lambda x: not self.item_limit or n <= self.item_limit, parser.parse(infilepath)):
                n += 1
                id = str(n)
                outdnadb.put(id, nseq['dna'])
                for (frame, dump, out) in translate_sequence(id, nseq['id'], nseq.get('desc', ''), nseq['dna']):
                    outprotdb.put(''.join([id, '_', str(frame)]), dump)
                    outfile.write(out)
            return n
        finally:
            outfile.close()

