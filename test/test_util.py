from util import *

def test_sequence_to_fasta():
    fasta = '>id123:abc_id id_id\nATTTACGPXYN\n'
    assert fasta == sequence_to_fasta('id123:abc_id id_id', 'ATTTACGPXYN')

def test_fix_fasta():
    assert '>id123:abc_id id_id\n' + 'A'*80 + '\n' + 'A'*80 + '\n' + 'A'*40 + '\n' == fixfasta('>id123:abc_id id_id\n' + 'A'*200 + '\n')
