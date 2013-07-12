from util.translator import *
import json

id = 'XYZ123'
name = 'id123:abc_id'
desc = 'kalaskul gen!!1'
seq = 'ATTTGCATTTTAGGGGCG'
revcomp = 'CGCCCCTAAAATGCAAAT'
prots = ['ICILGA', 'FAF*GR', 'LHFRGX', 'RP*NAN', 'APKMQX', 'PLKCKX']


def test_frame_sequence():
    for frame in range(0,3):
        assert seq[frame::] == frame_sequence(seq, frame)
        assert revcomp[frame::] == frame_sequence(seq, frame+3)

def test_translate_sequence():
    translation = translate_sequence(id, name, desc, seq)
    for frame in range(0,6):
        (tframe, js, fasta) = translation[frame]
        assert tframe == frame
        obj = json.loads(js)
        assert obj['protein'] == prots[frame]
        assert obj['name'] == name
        assert obj['description'] == desc
        assert obj['frame'] == frame



