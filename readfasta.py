#!/usr/bin/python
from __future__ import print_function
import time

import berkeley
from translator import translate_sequence

_DEBUG = True
_ITEM_LIMIT = 0

#Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
def translate_fasta(inpath, outpath):
    outfile = open(outpath, 'w')
    infile = open(inpath,'r')
    #outdb = berkeley.open_fragments('n')
    if _DEBUG:
        outdb = {}
        print('Start', time.asctime(time.localtime()))
    try:
        n = 0
        tempseq = []
        for line in infile:
            if line.startswith('>'):
                n += 1
                if _ITEM_LIMIT and n > _ITEM_LIMIT:
                    break
                if len(tempseq) > 0:
                    for (id, dump, out) in translate_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
                        outdb[id] = dump
                        outfile.write(out)
                (seqid, seqdesc) = line[1::].split(' ', 1)
                tempseq = []
            else:
                tempseq.append(line.rstrip())
        #When the file is finished: Save the final sequence just like the others
        for (id, dump, out) in translate_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
            outdb[id] = dump
            outfile.write(out)
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
    finally:
        infile.close()
        outfile.close()
        if _DEBUG:
            print('Finish:', time.asctime(time.localtime()))
        else:
            outdb.close()


if _DEBUG:
    import cProfile
    cProfile.run("translate_fasta('tutorial/database/ntsmall_plus_qnr.nfa', 'test.pfa')")
    #translate_fasta('tutorial/database/ntsubset_plus_7_qnr.nfa', 'test.pfa')
