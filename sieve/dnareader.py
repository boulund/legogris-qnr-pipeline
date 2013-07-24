#!/usr/bin/python
from __future__ import print_function
import types
import sys
from datetime import datetime
from itertools import takewhile

from parsing.fasta import FASTAParser
from parsing.fastq import FASTQParser
from util.translator import translate_sequence
from sieve import Sieve

class DNAReader(Sieve):
    """
    Reads in raw nucleotide sequences and translates them to protein sequences.
    FASTA and FASTQ input formats are supported in raw or gzip-compressed files.

    FASTA input::

        >id desc
        DNA

    DNA Database output::

        'nid' : 'DNA'

    Protein Database output::

        'nid_X' : JSON Serialization of translation output from `translator.translate_sequence` for frame X (0 <= X <= 5).

    Translated protein sequences for all six frames are written to output file in FASTA format.
    """
    def __init__(self, params, logfile):
        """
        Optional parameters:
            * item_limit (int, default 0): Maximum number of sequences to parse. 0 means no limit.
        """
        Sieve.__init__(self, params, logfile, name='FASTA translator', param_names=[('item_limit', 0)])

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

sieve = DNAReader
