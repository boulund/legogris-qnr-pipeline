from datetime import date
import time
import pickle

import fluff

TMPDIR = "./pipeline_data/"

class BLASTClusterer:
    def __init__(self, logfile):
        self.logfile = logfile

    # Number of CPUs to use, 0 means all
    # Percent identity threshold, range 3-100
    # Coverage threshold for blastclust, range 0.1-0.99
    def run(self, filepath, numcores, percent_identity, cov_threshold):
        logfile = self.logfile
        # Shorten the sequences that are too long, since blastclust seems to have
        # trouble aligning very long sequences (e.g. complete sequence genomes).
        try:
            SeqFilename = fluff.limit_sequence_length(filepath,64) # limit to 64 columns of sequence
        except fluff.PathError, e:
            logfile.write(e.message+"\n")
            logfile.write("\nThe clustering part of the pipeline is dependent of files from previous parts in the "+TMPDIR+" directory\n")
            exit(1)


        # Uniqueify the sequence IDs; needed for blastclust to cluster them since only
        # the first part of the identifier is parsed and thus sequence IDs risk becoming non-
        # unique. It is also a safeguard against redundant data sets.
        # (TODO: might be possible to use python "set" to remove duplicates, but
        # that would require checking of entire sequences to ensure robustness)
        # SeqFilename = TMPDIR+"retrieved_sequences.fasta.shortened"
        # UnSeqFilename = TMPDIR+"unique_retrieved_sequences.fasta.shortened"
        UnSeqFilename = SeqFilename+'.shortened'
        try:
            unique_sequences = fluff.uniqueify_seqids(SeqFilename,UnSeqFilename)
        except ValueError:
            logfile.write("Could not uniqueify the sequence IDs!\n")
            fluff.cleanup(TMPDIR)
            exit(1)


        # Run formatdb on the file outputted from uniqueify_seqids and
        # then run blastclust to cluster results (all in one function)
        t = time.asctime(time.localtime())
        logfile.write("Creating temporary database and running blastclust at: "+t+"\n")
        logfile.flush()

        try:
            #logfile.write("Running blastclust with the following settings:\n" \
            #              " Percent identity: "+str(percent_identity)+"\n Coverage Threshold: "+ \
            #              str(cov_threshold)+"\n Number of CPUs: "+str(numcores)+"\n")

            blastclust_return_text, blastclust_output = fluff.run_blastclust(UnSeqFilename,percent_identity,cov_threshold,numcores)

            if "error" in blastclust_return_text:
                logfile.write(blastclust_output+"\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            elif "ERROR" in blastclust_output[1]:
                logfile.write(blastclust_output[1]+"\n")
                fluff.cleanup(TMPDIR)
                exit(1)
            else:
                logfile.write(blastclust_return_text+"\n"+blastclust_output[0])

        except fluff.PathError, e:
            logfile.write(e.message+"\n"+UnSeqFilename+"\n")
            fluff.cleanup(TMPDIR)
            exit(1)


        # Parse output from blastclust (and de-uniqueify sequences IDs -- not needed anymore)
        blastclustoutputfile = UnSeqFilename+".clusters"
        try:
            parsedblastclust = fluff.parse_blastclust(blastclustoutputfile)
        except fluff.PathError, e:
            logfile.write(e.message+"\n"+blastclustoutputfile+"\n")
            fluff.cleanup(TMPDIR)
            exit(1)
        except ValueError:
            logfile.write("ERROR: Found nothing in blastclust output: "+blastclustoutputfile+"\n")
            fluff.cleanup(TMPDIR)
            exit(1)

        # Deunique:ify sequence IDs parsed from blastclust output,
        # needed only for writing out cluster scores later on
        clusters = fluff.deuniqueify_seqids(parsedblastclust)

        # Unpickle the scores_ids (needed for being able to run clustering separately)
        # Really unnecessary to do all the time but does not really matter
        try:
            pkfile = open(''.join([TMPDIR,"pickled.hsseq"]),'rb')
            scores_ids = pickle.load(pkfile)
        except IOError:
            logfile.write("Could not read the pickled high-scoring sequences (pickled.hsseq)\n")
            fluff.cleanup(TMPDIR)
        except pickle.UnpicklingError:
            logfile.write("Could not unpickle pickled.hsseq!\n")
            fluff.cleanup(TMPDIR)
        return (clusters, parsedblastclust, scores_ids)
