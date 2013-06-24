#!/usr/bin/python
from __future__ import print_function
import time
import uuid
import json
import random
from collections import defaultdict

import berkeley

_DEBUG = True
_ITEM_LIMIT = 10
_COMPLEMENTS = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'Y': 'R',
    'R': 'Y',
    'S': 'S',
    'W': 'W',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'V': 'B',
    'D': 'H',
    'H': 'D',
    'N': 'N',
    'X': 'X',
    '\n': ''
}
_GENCODE = defaultdict(lambda: 'X', {
    'TTT': 'F','TTC': 'F','TTY': 'F',
    'TTA': 'L','TTG': 'L','TTR': 'L',
    'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L','CTM': 'L','CTR': 'L','CTW': 'L','CTS': 'L','CTY': 'L','CTK': 'L','CTV': 'L','CTH': 'L','CTD': 'L','CTB': 'L','CTX': 'L','CTN': 'L',
    'ATT': 'I','ATC': 'I','ATA': 'I','ATY': 'I','ATW': 'I','ATM': 'I','ATH': 'I', 'ATG': 'M',
    'GTT': 'V','GTC': 'V','GTA': 'V','GTG': 'V','GTM': 'V','GTR': 'V','GTW': 'V','GTS': 'V','GTY': 'V','GTK': 'V','GTV': 'V','GTH': 'V','GTD': 'V','GTB': 'V','GTX': 'V','GTN': 'V',
    'TCT': 'S','TCC': 'S','TCA': 'S','TCG': 'S','TCM': 'S','TCR': 'S','TCW': 'S','TCS': 'S','TCY': 'S','TCK': 'S','TCV': 'S','TCH': 'S','TCD': 'S','TCB': 'S','TCX': 'S','TCN': 'S',
    'CCT': 'P','CCC': 'P','CCA': 'P','CCG': 'P','CCM': 'P','CCR': 'P','CCW': 'P','CCS': 'P','CCY': 'P','CCK': 'P','CCV': 'P','CCH': 'P','CCD': 'P','CCB': 'P','CCX': 'P','CCN': 'P',
    'ACT': 'T','ACC': 'T','ACA': 'T','ACG': 'T','ACM': 'T','ACR': 'T','ACW': 'T','ACS': 'T','ACY': 'T','ACK': 'T','ACV': 'T','ACH': 'T','ACD': 'T','ACB': 'T','ACX': 'T','ACN': 'T',
    'GCT': 'A','GCC': 'A','GCA': 'A','GCG': 'A','GCM': 'A','GCR': 'A','GCW': 'A','GCS': 'A','GCY': 'A','GCK': 'A','GCV': 'A','GCH': 'A','GCD': 'A','GCB': 'A','GCX': 'A','GCN': 'A',
    'TAT': 'Y','TAC': 'Y','TAY': 'Y',
    'TAA': '*','TAG': '*','TAR': '*',
    'CAT': 'H','CAC': 'H','CAY': 'H',
    'CAA': 'Q','CAG': 'Q','CAR': 'Q',
    'AAT': 'N','AAC': 'N','AAY': 'N',
    'AAA': 'K','AAG': 'K','AAR': 'K',
    'GAT': 'D','GAC': 'D','GAY': 'D',
    'GAA': 'E','GAG': 'E','GAR': 'E',
    'TGT': 'C','TGC': 'C','TGY': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R','CGC': 'R','CGA': 'R','CGG': 'R','CGM': 'R','CGR': 'R','CGW': 'R','CGS': 'R','CGY': 'R','CGK': 'R','CGV': 'R','CGH': 'R','CGD': 'R','CGB': 'R','CGX': 'R','CGN': 'R',
    'AGT': 'S','AGC': 'S','AGY': 'S',
    'AGA': 'R','AGG': 'R','AGR': 'R',
    'GGT': 'G','GGC': 'G','GGA': 'G','GGG': 'G','GGM': 'G','GGR': 'G','GGW': 'G','GGS': 'G','GGY': 'G','GGK': 'G','GGV': 'G','GGH': 'G','GGD': 'G','GGB': 'G','GGX': 'G','GGN': 'G',
    'CT': 'L', 'GT': 'V', 'TC': 'S', 'CC': 'P', 'AC': 'T', 'GC': 'A', 'CG': 'R', 'GG': 'G'
})

#Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
def translate_fasta(inpath, outpath):
    outfile = open(outpath, 'w')
    infile = open(inpath,'r')
    if _DEBUG:
        outdb = {}
        print('Start', time.asctime(time.localtime()))
    outdb = berkeley.open_fragments('n')
    try:
        n = 0
        tempseq = []
        for line in infile:
            if line.startswith('>'):
                n += 1
                if _ITEM_LIMIT and n > _ITEM_LIMIT:
                    break
                if len(tempseq) > 0:
                    for (id, dump, out) in _save_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
                        outdb[id] = dump
                        outfile.write(out)
                (seqid, seqdesc) = line[1::].split(' ', 1)
                tempseq = []
            else:
                tempseq.append(line.rstrip())
        #When the file is finished: Save the final sequence just like the others
        for (id, dump, out) in _save_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq)):
            outdb[id] = dump
            outfile.write(out)
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
    finally:
        infile.close()
        outfile.close()
        if _DEBUG:
            print('Finish:', time.asctime(time.localtime()))
        outdb.close()

#Translates the supplied DNA string in all 6 reading frames and stores the result in a FASTA format text file as well as in serialized JSON in a supplied key/value store.
def _save_sequence(name, desc, sequence):
    #Local variables = less overhead
    result = []
    gencode = _GENCODE
    complements = _COMPLEMENTS
    for frame in range(0,6):
        #First 3 frames are normal, following 3 are reverse complements
        if frame > 2:
            #Reverse and frame adjust
            dna = sequence[-(frame-2)::-1]
            #Complement
            dna  = ''.join([complements[c] for c in dna])
        else:
            dna  = sequence[frame::]
        #Translate dna codons to protein
        protein = ''.join([gencode[dna[i:i+3]] for i in xrange(0, len(dna), 3)])
        #Faster but less secure (wrt collissions) than stock uuid4
        id = uuid.UUID(int=random.getrandbits(128), version=4).hex
        seq = {
            'id': id,
            'dna': dna,
            'protein': protein,
            'name': name,
            'description': desc,
            'frame': frame+1
        }
        out = ''.join(['>', id, '\n', protein, '\n'])
        result.append((id, json.dumps(seq), out))
    return result

if _DEBUG:
    #import cProfile
    #cProfile.run("translate_fasta('tutorial/database/ntsmall_plus_qnr.nfa', 'test.pfa')")
    translate_fasta('tutorial/database/ntsubset_plus_7_qnr.nfa', 'test.pfa')
