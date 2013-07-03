from datetime import date
import time
import pickle
import json
from bsddb3 import db
import uuid

from sieve import Sieve
from fluff import PathError, cleanup
import fluff

def create(params, logfile):
    return BLASTClusterer(params, logfile)

class BLASTClusterer(Sieve):
    def init(self, params):
        self.indbflags = None
        self.indbmode = db.DB_HASH
        self.outdbflags = db.DB_DUPSORT
        self.outdbmode = db.DB_HASH
        self.name = 'BLASTClust'
        self.param_names = [
            'blastclust_out',
            ('clusters_out_path', ''),
            ('clusters_with_scores_out_path', ''),
            ('numcpu', 0), #Number of CPUs to use, 0 means all
            ('percent_identity', 90), # Percent identity threshold, range 3-100, default is 90 %
            ('coverage_threshold', 0.25), # Coverage threshold for blastclust, range 0.1-0.99
            ('reference_sequence_path', '') #Reference sequences to add to infile before clustering
        ]


    #def run(self, filepath, numcores, percent_identity, cov_threshold):
    #TODO: Ensure blastclust gets maximum 64 lines of FASTA (64*80 columns)
    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        logfile = self.logfile

        # Run formatdb on the file outputted from uniqueify_seqids and
        # then run blastclust to cluster results (all in one function)
        t = time.asctime(time.localtime())
        logfile.writeline("Running blastclust at: "+t)
        logfile.flush()

        if self.reference_sequence_path and path.isfile(self.reference_sequence_path):
            system('cat %s >> %s' % (self.reference_sequence_path, infilepath))
            logfile.write("Added reference sequences from "+self.reference_sequence_path+" to set to cluster\n")
        try:
            blastclust_return_text, blastclust_output, blastclustoutputfile = self.run_blastclust(infilepath)

            if "error" in blastclust_return_text:
                logfile.writeline(blastclust_output)
                exit(1)
            elif "ERROR" in blastclust_output[1]:
                logfile.writeline(blastclust_output[1])
                exit(1)
            else:
                logfile.writeline(blastclust_return_text)
                logfile.writeline(blastclust_output[0])

        except PathError, e:
            logfile.writeline(e.message)
            logfile.writeline(filename)
            exit(1)


        # Parse output from blastclust (and de-uniqueify sequences IDs -- not needed anymore)
        try:
            clusters = self.parse_blastclust(blastclustoutputfile)
        except PathError, e:
            logfile.write(e.message+"\n"+blastclustoutputfile+"\n")
            cleanup(TMPDIR)
            exit(1)
        except ValueError:
            logfile.write("ERROR: Found nothing in blastclust output: "+blastclustoutputfile+"\n")
            cleanup(TMPDIR)
            exit(1)

        # Output the identified cluster to files,
        # one with clean clusters and one with scores
        clusterout = withscores = None
        if self.clusters_out_path:
            clusterout = open(self.clusters_out_path,"w")
        if self.clusters_with_scores_out_path:
            withscores = open(self.clusters_with_scores_out_path,"w")
        outprotdb.truncate()
        try:
            for cluster in clusters:
                cid = uuid.uuid4().hex
                for seqID in cluster:
                    outprotdb.put(cid, seqID)
                    if clusterout or withscores:
                        seq = json.loads(inprotdb[seqID])
                        if clusterout:
                            clusterout.write(''.join([seq['name'],' ']))
                        if withscores:
                            withscores.write(''.join([seq['name'],"--",
                                        str(seq['score']),"--",str(seq['dscore']),' ']))
                if clusterout:
                    clusterout.write('\n')
                if withscores:
                    withscores.write('\n')
        finally:
            if clusterout:
                clusterout.close()
            if withscores:
                withscores.close()

        t = time.asctime(time.localtime())
        logfile.writeline("Found " + str(len(clusters)) + " clusters")
        logfile.writeline("Finished clustering sequences at: "+t)
        logfile.line()
        logfile.flush()
        return clusters

##-----------------------------------------------##
##                 RUN BLASTCLUST                ##
##-----------------------------------------------##
    def run_blastclust(self, infilepath):
        '''
        Run formatdb to create a BLAST database, then run
        blastclust on that database to cluster all hits.

        Input::

            infilepath  filename with sequences with unique ids to cluster.

        Returns::

            (None)      Writes output to a file, 'filename.clusters' that contains
                        all identified clusters on each row.

        Errors::

            PathError   raied if there is something wrong with the paths to output
                        or input files.
            ValueError  raised if there is something wrong

        '''

        from os import path
        import shlex, subprocess

        logfile = self.logfile
        if not(path.exists(path.abspath(infilepath))):
            raise PathError("ERROR: The path to unique sequence id hit file is incorrect")

        # Create string for calling the subprocess 'formatdb'
        formatdb_call = 'formatdb -t SignificantHits_UniqueSeqIDs -i %s' % infilepath
        try:
            formatdb = shlex.split(formatdb_call)
            subprocess.call(formatdb)
            return_text = "Created BLAST database with unique sequence IDs"
        except OSError:
            return_text = "The formatdb command could not be run:", formatdb_call

        # Run blastclust in a similar way
        outfilepath = ''.join([infilepath, '.clusters'])
        blastclust_call = 'blastclust -d %s -S %d -a %d -L %f -o %s' % (infilepath, self.percent_identity, self.numcpu, self.coverage_threshold, outfilepath)

        try:
            blastclust = shlex.split(blastclust_call)
            blastclust_output = subprocess.Popen(blastclust,\
                                    stdout=subprocess.PIPE,\
                                    stderr=subprocess.PIPE).communicate()
        except OSError:
            returnstring = ''.join(["ERROR: Blastclust could not be run. Is it properly installed?\n",
                                    "The following call was made: $", blastclust_call])
            return ("error",returnstring, outfilepath)


        return (return_text, blastclust_output, outfilepath)
############## END run_blastclust

##-----------------------------------------------##
##          PARSE OUTPUT FROM BLASTCLUST         ##
##-----------------------------------------------##
    def parse_blastclust(self, filename):
        '''
        Parses blastclust output into a nested list structure

        Input::

            filename    filename of blastclust output

        Returns::

            sequenceIDs list of sequence IDs with unique identifiers right after
                        the '>' symbol

        Errors::

            PathError   raised if the file does not exists
            ValueError  rasied if no unique identifiers could be found and removed

        '''

        from os import path
        logfile = self.logfile

        if not(path.isfile(path.abspath(filename))):
            raise PathError("ERROR: The blastclust output file could not be opened")

        try:
            filepath = path.abspath(filename)
            file = open(filepath,'r')
        except OSError:
            print "ERROR: Could not open", filepath
            exit(1)

        sequenceIDs = []

        try:
            for line in file:
                sequenceIDs.append(line.rstrip('\n ').split(' '))
        finally:
            file.close()

        if sequenceIDs == []:
            raise ValueError

        return sequenceIDs
############## END parse_blastclust

##-----------------------------------------------##
##              LIMIT SEQUENCE LENGTH            ##
##-----------------------------------------------##
def _limit_sequence_length(sequencefile,MAX_LINES=64):
    '''
    Takes a sequence file and shortens all sequences
    it to the MAX_LINES supplied.

    It is divided by
    inserting the sequence ID header ">seqid...."
    after the maximum sequence length position, thus
    creating several smaller sequence segments with
    the same sequence ID and header line.
    Recommended MAX_LINES is 64 (64 lines of 80
    amino acid residues = 5120).
    As the function splits at "\n" characters it
    requires the fasta file to keep sequences on
    several lines - fasta files with sequences on
    a single line will NOT work.

    Input::

        sequencefile    filename string.
        MAX_LINES       maximum number of lines for one sequence, defaults
                        to 64.

    Returns::

        filename (String)   Writes directly to disk to 'sequencefile.shortened'. filename is the filename of the output file.

    Errors::

        PathError   raised if there is something wrong with path

    '''

    from os import path

    if not path.isfile(path.abspath(sequencefile)):
        raise PathError(''.join(["ERROR: Path ",sequencefile," is not a valid file!"]))
    outfilename = sequencefile+'.shortened'

    try:
        file = open(path.abspath(sequencefile),'r')
        outfile = open(path.abspath(outfilename),'w')
    except OSError:
        raise PathError("Could not open sequence and/or output file!")

    # Reads the ENTIRE file into memory - it is assumed
    # that the number of sequences retrieved does not
    # require more memory than what can fit into memory.
    entire_file = file.read()

    # Splits at all sequences but the first (it does not start with a newline)
    # This method introduces an empty line between the first and following
    # sequences if it is split up, assumed non dangerous but still noted.
    sequences = entire_file.split("\n>")
    first = True # To give special treatment to first sequence
    seqout = []
    for sequence in sequences:
        if not sequence == "":
            if first:
                # Already has initial ">"
                lines = sequence.splitlines(True)
                seqout = []
                counter = 0
            else:
                # Add the'\n>' that was removed
                sequence = ''.join(["\n>",sequence])
                lines = sequence.splitlines(True)
                seqout = []
                counter = 0 # Reset, starting with new sequence

            for line in lines:
                # Add lines of sequence until MAX_LINES,
                # then add sequence identifier again and
                # continue until end of sequence
                if counter > MAX_LINES:
                    if first:
                        if not(lines[0].endswith("\n")):
                            seqout.append(''.join([lines[0],"\n"]))
                        else:
                            seqout.append(lines[0])
                    else:
                        # The following sequences
                        # got split at the initial '\n'
                        # an the identifier is thus in
                        # element 1 instead of 0
                        if not(lines[1].endswith("\n")):
                            seqout.append(''.join([lines[1],"\n"]))
                        else:
                            seqout.append(lines[1])
                    if not(line.endswith("\n")):
                        seqout.append(''.join([line,"\n"]))
                    else:
                        seqout.append(line)
                    counter = 0
                else:
                    seqout.append(line)
                    counter = counter + 1
            #seqout.append("") # Probably not needed...
            first = False # Now the first has been processed

            outfile.write(''.join(seqout))

    return outfilename
############## END limit_sequence_length

