#!/usr/bin/python
from os import path
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import uuid
import json

import berkeley

def translate_fasta(inpath, outpath):
    outfile = open(outpath, 'w')
    db = berkeley.open()
    try:
        for r in SeqIO.parse(inpath, 'fasta'):
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