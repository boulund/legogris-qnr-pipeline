#!/usr/bin/python
from __future__ import print_function
from os import path
import os
import select
import time
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import uuid
import json
import subprocess, shlex
import re

import berkeley

_ITEM_LIMIT = 0
_COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#TODO: codontablegen.py
_GENCODE = {
  re.compile('TT[TCY]'): 'F',
  re.compile('TT[AGR]'): 'L',
  re.compile('CT.'): 'L',
  re.compile('AT[TCAYWMH]'): 'I',
  re.compile('ATG'): 'M',
  re.compile('GT.'): 'V',
  re.compile('TC.'): 'S',
  re.compile('CC.'): 'P',
  re.compile('AC.'): 'T',
  re.compile('GC.'): 'A',
  re.compile('TA[TCY]'): 'Y',
  re.compile('TA[AGR]'): '*',
  re.compile('CA[TCY]'): 'H',
  re.compile('CA[AGR]'): 'Q',
  re.compile('AA[TCY]'): 'N',
  re.compile('AA[AGR]'): 'K',
  re.compile('GA[TCY]'): 'D',
  re.compile('GA[AGR]'): 'E',
  re.compile('TG[TCY]'): 'C',
  re.compile('TGA'): '*',
  re.compile('TGG'): 'W',
  re.compile('CG.'): 'R',
  re.compile('AG[TCY]'): 'S',
  re.compile('AG[AGR]'): 'R',
  re.compile('GG.'): 'G'
}

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
    outdb = berkeley.open_fragments('n')
    db = {} #berkeley.open_dna_input('n')
    outdb = berkeley.open_fragments('n')
    #p = subprocess.Popen(args, stdin=subprocess.PIPE,stdout=subprocess.PIPE,bufsize=-1)
    #lr = LineReader(p.stdout, p.stdin, db, outdb)
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
                    for frame in range(1,4):
                        dna = _construct_sequence(seqid, ''.join(tempseq), frame, False)
                        protein = _translate_dna(dna)
                        seq = {
                            'id': uuid.uuid4().hex,
                            'dna': dna,
                            'protein': protein,
                            'name': seqid,
                            'description': seqdesc.lstrip(),
                            'frame': int(frame)
                        }
                        outdb[seqid] = json.dumps(seq)

                (seqid, seqdesc) = line[1::].split(' ', 1)
                tempseq = []
            else:
                tempseq.append(line.rstrip())
        #TODO: Save the last one too
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
    finally:
        infile.close()
        outfile.close()
        outdb.close()

fasta_regex = re.compile('>([^\s]*_([0-6]))([^\r\n]*)\n([^>]+)')

def _save_sequence(name, desc, dna, db, transeq):
    s = ''.join(['>',name,' ', desc, '\n',dna,'\n'])
    ids = name.rsplit('|')
    ids.reverse()
    id = ''
    for _id in ids:
        id = _id #transeq extracts id from id string, so do the same thing
        if id != '':
            break
    print('Saving ', id, ' from ', name)
    db[id] = dna
    #transeq.stdin.write(s)
    #transeq.stdin.flush()
    #result = p.communicate(input=s)[0]
    #result = p.stdout.read(-1)
'''    for (name, frame, desc, protein) in fasta_regex.findall(result):
        id = uuid.uuid4().hex
        seq = {
            'id': id,
            'dna': dna,
            'protein': protein,
            'name': name,
            'description': desc.lstrip(),
            'frame': int(frame)
        }
        #db[seq['id']] = json.dumps(seq)
        out = ''.join(['>', id, '\n', protein, '\n'])
        outfile.write(out) '''

def _construct_sequence(name, dna, frame=1, reverse=False):
    if reverse:
        dna = _reverse_complement(dna)[frame-1::]
        frame += 3
    else:
        dna = dna[frame-1::]
    return dna

def _translate_dna(dna):
    pseq = []
    for i in xrange(0, len(dna), 3):
        codon = dna[i:i+3]
        match = False
        for (rex, acid) in _GENCODE.items():
            if rex.match(codon):
                pseq.append(acid)
                match = True
                break
        if not match:
            pseq.append('X')
    return ''.join(pseq)

def _reverse_complement(dna):
    cs = []
    for c in reversed(dna):
        if c != '\n':
            cs.append(_COMPLEMENTS[c])
    return ''.join(cs)

translate_fasta('tutorial/database/ntsmall_plus_qnr.nfa', 'test.pfa')
