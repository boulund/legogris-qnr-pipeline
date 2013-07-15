from __future__ import print_function
import shlex
import subprocess
import json

from parsing.fasta import FASTAParser
from util import sequence_to_fasta
from sieve import Sieve
from util import PathError

def create(params, logfile):
    return SGAAligner(params, logfile)

class SGAAligner(Sieve):
    def init(self, params):
        self.indbmode = True
        self.outdbmode = True
        self.name = 'SGA Aligner'
        self.param_names = [
            ('result_out_path', ''),
            ('max_edges', 200),
            ('numcpu', 4),
            ('min_branch_length', 45),
            ('min_merge_overlap', 15),
            ('min_assembly_overlap', 20),
            ('max_gap_divergence', 0),
            ('resolve_small', 10),
            ('error_rate', 0.02)
        ]


    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        # Run SGA pipeline
        try:
            (dir, filepath) = infilepath.rsplit('/', 1)
            commands = [
                str.format('sga preprocess --pe-mode=0 -o {0}.processed --permute-ambiguous {0}', infilepath),
                str.format('sga index -a ropebwt -t {0} --no-reverse {1}.processed', self.numcpu, infilepath),
                str.format('sga correct -a overlap -t {0} -o {1}.corrected {1}.processed', self.numcpu, infilepath),
                str.format('sga index -a ropebwt -t {0} {1}.corrected', self.numcpu, infilepath),
                str.format('sga filter -t {0} --no-kmer-check {1}.corrected', self.numcpu, infilepath),
                str.format('sga fm-merge -m {0} -t {1} -o {2}.merged {2}.filter.pass.fa', self.min_merge_overlap, self.numcpu, infilepath),
                str.format('sga index -d 1000000 -t {0} {1}.merged', self.numcpu, infilepath),
                str.format('sga rmdup -t {0} -o {1}.merged.rmdup.fa {1}.merged', self.numcpu, infilepath),
                str.format('sga overlap --min-overlap={0} -e {1} --exhaustive --threads={2} {3}.merged.rmdup.fa', self.min_merge_overlap, self.error_rate, self.numcpu, infilepath),
                str.format('assemble --max-edges={0} --min-branch-length={1} --min-overlap={2} --max-gap-divergence={3} --resolve-small={4} -o {5}.result {5}.merged.rmdup.asqg.gz', self.max_edges, self.min_branch_length, self.min_assembly_overlap, self.error_rate, self.resolve_small, self.numcpu, infilepath)
            ]
            for cmd in commands:
                call = shlex.split(cmd)
                output = subprocess.Popen(call, stdin=subprocess.PIPE,
                                            stderr=subprocess.PIPE).communicate()
                print(output[1])
            self.logfile.write("Finished SGA alignment on file "+infilepath+"\n")
        except OSError, e:
            self.logfile.write("Could not open one of sga align or "+infilepath+"\n")
            raise e
        self.logfile.flush()
        if self.parse_output:
            #There are two sets of output files: One with the final aligned contigs and one with the removed duplicates.
            #Step one: Associate aligned contigs with original fragments
            parser = FASTAParser(self.logfile)
            aligned = [seq for seq in parser.parse_fasta(infilepath+'.result-contigs.fa')]
            outfile = open(outfilepath, 'w')
            try:
                for pseq in parser.parse_fasta(infilepath+'.filter.pass.fa'):
                    seq = json.loads(inprotdb.get(pseq['id']))
                    for aseq in aligned:
                        if pseq['dna'] in aseq['dna']:
                            seq['contig'] = aseq['id']
                            outprotdb.put(pseq['id'], json.dumps(seq))
                            outfile.write(sequence_to_fasta(''.join([pseq['id'], '_', seq['name'], '_', seq['contig']]), seq['protein']))
            finally:
                outfile.close()
