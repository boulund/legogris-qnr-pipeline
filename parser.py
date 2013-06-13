from os import path, makedirs, system
from datetime import date
import time
import pickle

import fluff

logfileseparator = "----------------------------------------------------------------------"
TMPDIR = "./pipeline_data/"

class Parser:
    def __init__(self, logfile, classificationfunction):
        self.classificationfunction = classificationfunction
        self.logfile = logfile

    def parse_files(self, files, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright):
        logfile = self.logfile
        t = time.asctime(time.localtime())
        print "Starting to extract and classify hits from hmmsearch output at: "+t
        logfile.write("Starting to extract and classify hits from hmmsearch output at: "+t+"\n")

        # Check that all paths given are valid and that files exists in those locations
        hmmsearch_result_files = []
        for filepath in files:
            if path.isfile(path.abspath(filepath)):
                hmmsearch_result_files.append(path.abspath(filepath))
            else:
                print " ERROR: incorrect path for hmmsearch output file:",filepath
                logfile.write(" ERROR: incorrect path for hmmsearch output file: "+filepath+"\n")

        # Open file for writing the retrieved sequences to semi-temporary file
        try:
            retrseqfile = open(outfile,'w')
        except OSError:
            print " ERROR: Could not open file for writing:", RETR_SEQ_FILEPATH
            logfile.write(" ERROR: Could not open file for writing: "+RETR_SEQ_FILEPATH+"\n")


        numerrors = 0
        scores_ids = []
        for hmmsearch_result_file in hmmsearch_result_files:
            sequences, score_id_tuples = self.parse_file(hmmsearch_result_file, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright)
            if sequences:
                for sequence in sequences:
                    retrseqfile.write(''.join([sequence,"\n"]))
            scores_ids.append(score_id_tuples)
            logfile.flush()

        if numerrors >= len(hmmsearch_result_files):
            print "CATASTROPHIC: No input files contains any potential hits!?"
            logfile.write("CATASTROPHIC: No input files contains any potential hits!?\n")
            exit(1)

        # Pickle the scores_ids for disconnected usage in other scripts
        # note that scores_ids is a very messy structure!
        data = scores_ids
        pickle_filename = ''.join([TMPDIR,"pickled.hsseq"])
        outpickle = open(pickle_filename,'wb')
        pickle.dump(data,outpickle)
        outpickle.close()
        retrseqfile.close()

        t = time.asctime(time.localtime())
        print "Finished parsing hmmsearch output files and classifying hits at: "+t
        logfile.write("Finished parsing hmmsearch output files and classifying hits at: "+t+"\n")
        print logfileseparator
        logfile.write(logfileseparator+"\n")
        logfile.flush()

    def parse_file(self, hmmsearch_result_file, outfile, minscore, retrdb, classifyC, classifyD, extendleft, extendright):
        logfile = self.logfile
        classificationfunction = self.classificationfunction
        # Open file for reading
        try:
            # Parse the hmmsearch output file -- Extract hits above minscore (0)
            print "Parsing",hmmsearch_result_file
            logfile.write("Parsing "+hmmsearch_result_file+"\n")
            parsed = fluff.parse_hmmsearch_output(hmmsearch_result_file,minscore)
            score_id_tuples, dbpath = parsed # Unpack parsed information
            scores,dscores,ids = zip(*score_id_tuples) # Unzip the scores/IDs
            errmessages = []
            if retrdb:
                # Retrieve hits directly from their source database, if
                # they classify correctly according to the classification
                # function. The database needs an index file for this,
                # created using cdbfasta (standard settings), if not available
                # things will go bad.
                print "Retrieving full length sequences from database"
                logfile.write("Retrieving full length sequences from database")
                sequences, errmessages = fluff.retrieve_sequences_from_db(dbpath, ids, dscores, outfile,
                                                func=classificationfunction,
                                                longseqcutoff=classifyC,
                                                longseqdef=classifyD)
            elif int(extendleft) or int(extendright):
                #print "EXTENDING HITS" # DEBUG
                # Extend the hits to the HMMER model and retrieve "a little more"
                # around the edges of the hits, adds the sequence origin to the
                # fasta headers. Only retrieves sequences that classify correctly
                # according to the classification function.
                print "Retrieving extended sequences from database that classified as potential hits"
                logfile.write("Retrieving extended sequences from database that classified as potential hits")
                sequences, errmessages = fluff.extend_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                ids, minscore, dbpath,
                                                                extendleft=extendleft,
                                                                extendright=extendright,
                                                                func=classificationfunction,
                                                                longseqcutoff=classifyC,
                                                                longseqdef=classifyD)
            else:
                # Retrieve interesting sequences from hmmsearch output
                # and write out to RETR_SEQ_FILEPATH
                # This is where the sequences get their origin attached
                # Classify the hmmsearch hits according to classification function
                # and only write sequences to disk if they are classified as 'true'
                print "Retrieving sequences from hmmsearch output that classified as potential hits"
                logfile.write("Retrieving sequences from hmmsearch output that classified as potential hits")
                sequences = fluff.retrieve_sequences_from_hmmsearch(hmmsearch_result_file,
                                                                    ids, minscore, dbpath,
                                                                    func=classificationfunction,
                                                                    longseqcutoff=classifyC,
                                                                    longseqdef=classifyD)

            if errmessages:
                for message in errmessages:
                    print message,
                    logfile.write(message)
            print "Retrieved "+str(len(sequences))+" full length sequences from database"
            logfile.write("Retrieved "+str(len(sequences))+" full length sequences from database\n")
            # Write the identified sequences (fragments/domains) to disk,
            # they have been classified inside the previous function and
            # non-qnr like sequences have been removed.
            return (sequences, score_id_tuples)

        except IOError:
            print "Could not open file for reading:", hmmsearch_result_file
            print "Continuing with next file..."
            logfile.write("Could not open file for reading: "+hmmsearch_result_file+"\n")
            logfile.write("Continuing with next file...\n")
        except ValueError:
            print "Found no sequences with domain score above or equal", minscore,
            print " in file:",hmmsearch_result_file
            logfile.write("Found no sequences with domain score above or equal "+str(minscore))
            logfile.write(" in file: "+hmmsearch_result_file+"\n")
        except fluff.PathError, e:
            print e.message
            logfile.write(e.message+"\n")
        except fluff.ParseError as e:
            print e.message
            logfile.write(e.message+"\n")
        return ([], score_id_tuples)

