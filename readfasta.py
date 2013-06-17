#!/usr/bin/python
from os import path
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import uuid
import json

import berkeley

ITEM_LIMIT = 0

#Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
def translate_fasta(inpath, outpath):
    outfile = open(outpath, 'w')
    db = berkeley.open_fragments('n')
    try:
        n = 0
        for r in SeqIO.parse(inpath, 'fasta'):
            n += 1
            if ITEM_LIMIT and n > ITEM_LIMIT:
                break
            for j in range(0,2):
                for i in range(0,3):
                    dna = r.seq[i::] if j else r.seq.reverse_complement()[i::]
                    seq = {'id': uuid.uuid4().hex,
                        'dna': str(dna),
                        'protein': str(dna.translate(table=11)),
                        'name': r.name,
                        'frame': i+1 + 3*j
                        }
                    db[seq['id']] = json.dumps(seq)
                    outfile.write('>')
                    outfile.write(seq['id'])
                    outfile.write("\n")
                    outfile.write(seq['protein'])
                    outfile.write("\n")
        return outpath
    finally:
        outfile.close()
        db.close()

##translate_fasta('tutorial/database/database.nfa', 'test.pfa')
