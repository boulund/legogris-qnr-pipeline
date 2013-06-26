#!/usr/bin/python
from __future__ import print_function
import types
import sys

from bsddb import db
from translator import translate_sequence

_DEBUG = False
_ITEM_LIMIT = 0

class Sieve(object):
    def __init__(self, params):
        self.outdbflags = None
        self.outdbmode = db.DB_HASH
        self.name = 'FASTA translator'
        self.params = params

    def run(self, indb, infile, outdb, outfile):
        n = 0
        tempseq = []
        for line in infile:
            if line.startswith('>'):
                n += 1
                if _ITEM_LIMIT and n > _ITEM_LIMIT:
                    break
                if len(tempseq) > 0:
                    for (id, dump, out) in translate_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
                        outdb.put(id, dump)
                        outfile.write(out)
                (seqid, seqdesc) = line[1::].split(' ', 1)
                tempseq = []
            else:
                tempseq.append(line.rstrip())
        #When the file is finished: Save the final sequence just like the others
        for (id, dump, out) in translate_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
            outdb.put(id, dump)
            outfile.write(out)


if __name__ == "__main__":
    import cProfile
    import berkeley
    #Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
    #TODO: Clean up
    def translate_fasta(inpath, outpath):
        outfile = open(outpath, 'w')
        infile = open(inpath,'r')
        outdb = berkeley.open_fragments('n')
        Sieve().run(None, infile, outdb, outfile)
    #inpath = 'tutorial/database/ntsmall_plus_qnr.nfa'
    inpath = 'tutorial/database/ntsubset_plus_7_qnr.nfa'
    cProfile.run("translate_fasta('"+inpath+"', 'test.pfa')")
    #translate_fasta(inpath, 'test.pfa')
