from __future__ import print_function
import shlex
import subprocess
import json
import os

from parsing.fasta import FASTAParser
from util import sequence_to_fasta
from sieve import Sieve
from util import PathError

class SGAAligner(Sieve):
    """
    Run the String Graph Assembler pipeline (preprocess, correct, filter, fm-merge, rmdup, overlap, assemble).
    Output the assembled contigs.
    """

    def __init__(self, params, logfile):
        """
        Optional parameters:
            * parse_output (bool, True): Whether the sieve should write its own output file.

        Optional parameters passed to SGA:
            * max_edges (int, 400)
            * numcpu (int, 4)
            * min_branch_length (int, 35)
            * min_merge_overlap (int, 15)
            * min_assembly_overlap (int, 5)
            * max_gap_divergence (int, 0)
            * resolve_small (int, 500)
            * error_rate (float, 0.02)
        """

        param_names = [
            ('parse_output', True),
            ('max_edges', 400),
            ('numcpu', 4),
            ('min_branch_length', 35),
            ('min_merge_overlap', 15),
            ('min_assembly_overlap', 5),
            ('max_gap_divergence', 0),
            ('resolve_small', 500),
            ('error_rate', 0.02)
        ]
        Sieve.__init__(self, params, logfile, name='SGA', param_names=param_names)


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        # Run SGA pipeline
        (dir, outfilename) = outfilepath.rsplit('/', 1)
        (indir, infilename) = infilepath.rsplit('/', 1)
        filepath = '/'.join([dir, infilepath.rsplit('/', 1)[1]])
        commands = [
            str.format('sga preprocess --pe-mode=0 -o {0}.processed --permute-ambiguous {1}', infilename, os.path.abspath(infilepath)),
            str.format('sga index -a ropebwt -t {0} --no-reverse {1}.processed', self.numcpu, infilename),
            str.format('sga correct -a overlap -t {0} -o {1}.corrected {1}.processed', self.numcpu, infilename),
            str.format('sga index -a ropebwt -t {0} {1}.corrected', self.numcpu, infilename),
            str.format('sga filter -t {0} --no-kmer-check --substring-only {1}.corrected', self.numcpu, infilename),
            str.format('sga fm-merge -m {0} -t {1} -o {2}.merged {2}.filter.pass.fa', self.min_merge_overlap, self.numcpu, infilename),
            str.format('sga index -d 1000000 -t {0} {1}.merged', self.numcpu, infilename),
            str.format('sga rmdup -t {0} -o {1}.merged.rmdup.fa {1}.merged', self.numcpu, infilename),
            str.format('sga overlap --min-overlap={0} -e {1} --exhaustive --threads={2} {3}.merged.rmdup.fa', self.min_merge_overlap, self.error_rate, self.numcpu, infilename),
            str.format('sga assemble --max-edges={0} --min-branch-length={1} --min-overlap={2} --max-gap-divergence={3} --resolve-small={4} -o {5}.result {5}.merged.rmdup.asqg.gz', self.max_edges, self.min_branch_length, self.min_assembly_overlap, self.max_gap_divergence, self.resolve_small, infilename)
        ]
        for cmd in commands:
            call = shlex.split(cmd)
            output = subprocess.Popen(call, cwd=dir, stdin=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            #print(output[1])
        self.logfile.write("Finished SGA alignment on file "+infilepath+"\n")
        self.logfile.flush()
        if self.parse_output:
            #There are two sets of output files: One with the final aligned contigs and one with the removed duplicates.
            #Step one: Associate aligned contigs with original fragments
            #This does obviously not work with mutation and error rate > 0. Should probably use MAFFTA or something like that instead
            parser = FASTAParser(self.logfile)
            aligned = [seq for seq in parser.parse_fasta(dir + '/' + infilename +'.result-contigs.fa')]
            outfile = open(outfilepath, 'w')
            try:
                for pseq in parser.parse_fasta(dir + '/' + infilename +'.filter.pass.fa'):
                    seq = json.loads(inprotdb.get(pseq['id']))
                    for aseq in aligned:
                        if pseq['dna'] in aseq['dna']:
                            seq['contig'] = aseq['id']
                            outprotdb.put(pseq['id'], json.dumps(seq))
                            #outfile.write(sequence_to_fasta(''.join([pseq['id'], '_', seq['name'], '_', seq['contig']]), seq['protein']))
            finally:
                outfile.close()

sieve = SGAAligner
