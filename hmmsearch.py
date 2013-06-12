from os import path, makedirs, system
import shlex
import subprocess
from datetime import date
import time
logfileseparator = "----------------------------------------------------------------------"
class HMMSearch:
    def __init__(self, logfile):
        self.logfile = logfile
    def search(
            self,
            model,# Path to MODEL
            hmmsearch_outdir,# Path to output directory
            numcpu,# CPU flag; "--cpu 4" means four CPUs, empty one core
            wrap_long_lines,
            use_heuristics,
            databases
        ):
        logfile = self.logfile
        cpuflag = ''.join(["--cpu ",str(numcpu)])
        # Heuristics on/off; --max means no heuristics (max sensitivity), empty full heuristics
        # There is little reason not to use heuristics, HMMer has a higher propensity
        # for crashing if not used and it only increases the number of really low
        heurflag = '--max' if use_heuristics else ''
        # Text wrap long lines; "--notextw" means true, empty false
        # THIS MUST BE ENABLED TO ENSURE CORRECT BEHAVIOR!
        textwflag = '--notextw' if wrap_long_lines else ''

        args = []

        # Retrieve current date, used in output filename to unique:ify the output filenames
        d = date.today()
        t = time.asctime(time.localtime())

        ## HMMsearch with the settings defined above
        # Print a log of the hmmsearch run settings
        logfile.write("Running hmmsearch at:"+t+"\n")
        print "Running hmmsearch at:",t
        #logfile.write("Model used                       : "+path.basename(model)+"\n")
        #logfile.write("Output directory                 : "+hmmsearch_outdir+"\n")
        #logfile.write("CPU-flag (no flag means one cpu) : "+cpuflag+"\n")
        #logfile.write("Text Wrap (empty means textwrap) : "+textwflag+"\n")
        #logfile.write("Sensitivity (empty means default): "+heurflag+"\n")
        #logfile.write("\nThe following input file(s) were entered at command line:\n"+'\n'.join(databases)+"\n")
        #print "Model used                       :", path.basename(model)
        #print "Output directory                 :", hmmsearch_outdir
        #print "CPU-flag (no flag means one cpu) :", cpuflag
        #print "Text Wrap (empty means textwrap) :", textwflag
        #print "Sensitivity (empty means default):", heurflag
        #print "\nThe following input file(s) were entered at command line:\n", '\n'.join(databases)

        for database in databases:
            # Retrieve the filename and path of the database
            outfilename = path.basename(database)
            database = path.abspath(database)

            # Put together the entire string to call hmmsearch
            call_list = ''.join(["hmmsearch ", cpuflag, " ", textwflag, " ", heurflag, " ",
                                "-o ", hmmsearch_outdir, outfilename,
                                ''.join([".hmmsearched--", d.isoformat(), " "]),
                                model, " ", database])
            hmmsearch = shlex.split(call_list)
            # Run hmmsearch
            try:
                output = subprocess.Popen(hmmsearch, stdin=subprocess.PIPE,
                                            stderr=subprocess.PIPE).communicate()
                if "Error: Failed to open hmm file" in output[1]:
                    print "CATASTROPHIC: Could not open HMM:",model
                    print "Make sure 'model.hmm' is available in current directory or"
                    print "supply the -m argument with path to your HMM file"
                    logfile.write("CATASTROPHIC: Could not open HMM: "+model+"\n")
                    logfile.write("Make sure 'model.hmm' is available in current directory or\n")
                    logfile.write("supply the -m argument with path to your HMM file\n")
                    exit(1)
                if use_heuristics:
                    args.append(hmmsearch[6]) # 5 contains the output file path, used later
                else:
                    args.append(hmmsearch[5]) # 5 contains the output file path, used later
                print "Finished hmmsearch on file",database
                logfile.write("Finished hmmsearch on file "+database+"\n")
            except OSError:
                print "Could not open:", database
                logfile.write("Could not open: "+database+"\n")
            logfile.flush()

        # Output some more details for the log
        t = time.asctime(time.localtime())
        print "Finished hmmsearching the databases at:", t
        print logfileseparator
        logfile.write("Finished hmmsearching the databases at: "+t+"\n")
        logfile.write(logfileseparator+"\n")
        logfile.flush()
        return args
if __name__ == 'main':
    pass
