#!/usr/bin/python
from __future__ import print_function
import types
import sys
import leveldb
from datetime import datetime

from translator import translate_sequence
from sieve import Sieve

def create(params, logfile):
    return FastaReader(params, logfile)

class FastaReader(Sieve):
    def init(self, params):
        self.outdbmode = True
        self.name = 'FASTA translator'
        self.param_names = [('item_limit', 0), ('db_batch_size', 100000)]

    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        infile = open(infilepath,'r')
        outfile = open(outfilepath, 'w')
        dnabatch = leveldb.WriteBatch()
        protbatch = leveldb.WriteBatch()
        try:
            n = 0
            tempseq = []
            for line in infile:
                if line.startswith('>'):
                    n += 1
                    if self.item_limit and n > self.item_limit:
                        break
                    if len(tempseq) > 0:
                        dna = ''.join(tempseq)
                        id = str(n)
                        seqs = translate_sequence(id, seqid, seqdesc.lstrip(), dna)
                        dnabatch.Put(id, dna)
                        for (frame, dump, out) in seqs:
                            protbatch.Put(''.join([id, '_', str(frame)]), dump)
                            outfile.write(out)
                        if n % self.db_batch_size == 0:
                            flushstarttime = datetime.now()
                            self.logfile.writeline("Flushing to disk at %s" % flushstarttime)
                            outdnadb.Write(dnabatch) #, sync=True
                            outprotdb.Write(protbatch) #, sync=True
                            dnabatch = leveldb.WriteBatch()
                            protbatch = leveldb.WriteBatch()
                            self.logfile.writeline("Done flushing to disk: %s" % (datetime.now() - flushstarttime))
                            self.logfile.writeline(outdnadb.GetStats())
                            self.logfile.writeline(outprotdb.GetStats())
                    (seqid, seqdesc) = line[1::].split(' ', 1)
                    tempseq = []
                else:
                    tempseq.append(line.rstrip())
            #When the file is finished: Save the final sequence just like the others
            dna = ''.join(tempseq)
            id = str(n)
            seqs = translate_sequence(id, seqid, seqdesc.lstrip(), dna)
            dnabatch.Put(id, dna)
            for (frame, dump, out) in seqs:
                protbatch.Put(''.join([id, '_', str(frame)]), dump)
                outfile.write(out)
        finally:
            flushstarttime = datetime.now()
            self.logfile.writeline("Flushing to disk at %s" % flushstarttime)
            outdnadb.Write(dnabatch) #, sync=True
            outprotdb.Write(protbatch) #, sync=True
            self.logfile.writeline("Done flushing to disk: %s" % (datetime.now() - flushstarttime))
            self.logfile.writeline(outdnadb.GetStats())
            self.logfile.writeline(outprotdb.GetStats())
            infile.close()
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
