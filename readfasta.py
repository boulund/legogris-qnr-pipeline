#!/usr/bin/python
from __future__ import print_function
from os import path
import os
import select
import time
import uuid
import json
import subprocess, shlex
import re
from collections import defaultdict

import berkeley

_ITEM_LIMIT = 0
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

class LineReader(object):
    def __init__(self, fd, transeq, indb, outdb):
        self.transeq = transeq
        self._fd = fd
        self._poll = select.poll()
        self._poll.register(fd)
        self.frames = []
        self.tempseq = []
        self.seqid = ''
        self.seqdesc = ''
        self.indb = indb
        self.outdb = outdb

    def process(self):
        while self._poll.poll(0):
            self._process_line(self._fd.readline())

    def _process_line(self, line):
        if line.startswith('>'):
            if len(self.tempseq) > 0:
                self.frames.append((self.seqid, self.seqdesc, ''.join(self.tempseq)))
                dna = ''
                if len(self.frames) == 6:
                    print('id', self.seqid)
                    name = self._save_translations()
                    self.frames = []
                    del self.indb[name]
            (self.seqid, self.seqdesc) = line[1:-1].split(' ', 1)
            self.tempseq = []
        else:
            self.tempseq.append(line.rstrip())

    def _save_translations(self):
        result = []
        dna = ''
        name = ''
        description = ''
        try:
            for (parse_id, description, protein) in self.frames:
                (name, frame) = parse_id.rsplit('_', 1)
                if dna == '':
                    dna = self.indb[name]
                seq = {
                    'id': uuid.uuid4().hex,
                    'dna': dna,
                    'protein': protein,
                    'name': name,
                    'description': description,
                    'frame': int(frame)
                }
        except KeyError:
            print(self.frames)
            print(self.indb)
        return name

#Reads FASTA file. Adds all six frames of fragments to database and a new FASTA file with new UUIDs as keys in both.
def translate_fasta(inpath, outpath):
    outfile = open(outpath, 'w')
    #outdb = berkeley.open_fragments('n')
    outdb = {} #Just for profiling/testing: No db
    infile = open(inpath,'r')
    #First step: Parse fasta file, save fragments to db and pipe to transeq process
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
                    _save_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq), outdb, outfile)
                (seqid, seqdesc) = line[1::].split(' ', 1)
                tempseq = []
            else:
                tempseq.append(line.rstrip())
        _save_sequence(seqid, seqdesc.lstrip(), ''.join(tempseq), outdb, outfile)
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
    finally:
        infile.close()
        outfile.close()
        #outdb.close()
        print('Finish:', time.asctime(time.localtime()))

def _save_sequence(name, desc, sequence, outdb, outfile):
    #Reverse complement
    revseq = ''.join([_COMPLEMENTS[c] for c in sequence[::-1]])
    for frame in range(0,6):
        if frame > 2:
            dna = revseq[frame-3::]
        else:
            dna = sequence[frame::]
        #translate dna codons to protein
        protein = ''.join([_GENCODE[dna[i:i+3]] for i in xrange(0, len(dna), 3)])
        id = uuid.uuid4().hex
        seq = {
            'id': id,
            'dna': dna,
            'protein': protein,
            'name': name,
            'description': desc,
            'frame': frame
        }
        outdb[id] = json.dumps(seq)
        out = ''.join(['>', id, '\n', protein, '\n'])
        outfile.write(out)

import cProfile

cProfile.run("translate_fasta('tutorial/database/ntsmall_plus_qnr.nfa', 'test.pfa')")
#translate_fasta('tutorial/database/ntsmall_plus_qnr.nfa', 'test.pfa')
