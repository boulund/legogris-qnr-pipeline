README
Qnr-search pipeline 
version 0.8055 beta

Fredrik Boulund
2012-05-02
fredrik.boulund@chalmers.se
Chalmers University of Technology

___________________________
TABLE OF CONTENTS
1.      Introduction
2.      Installation
2.1     Required packages
3.      The pipeline
3.1     Pipeline function
3.1.1   Pipeline structure/layout
3.1.2   Command line options
3.1.3   Pipeline output folder structure
3.2     Database preparation
3.3     Usage example
4.      Known issues
5.      Source code
6.      Acknowledgements
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯


1. Introduction
This pipeline was developed during the summer of 2010 as a project in the later
stages of my Master of Science degree in Bioinformatics and Systems Biology at 
Chalmers University, under the supervision of Erik Kristiansson and Anna 
Johnning, both at Sahlgrenska Academy. Development was continued at Chalmers 
University of Technology in 2011. The original purpose of the pipeline was to 
search large amounts of metagenomic data using the latest available version of 
HMMER to leverage the new speed improvements of this package and the 
specificity/sensitivity of using a Hidden Markov Model in the search for novel 
antibiotic resistance genes.



2. Installation
Note that as the pipeline relies on several components of a standard Linux 
distribution it only works on Linux systems, preferably with as much RAM and 
hard drive space available as possible. It can easily consume RAM in excess of
12GB and will require at least 3-4 times as much storage space as the 
untranslated nucleotide database you are interested in searching. Being a 
Python program, the pipeline naturally requires Python to run (it should be 
included in any standard linux distribution). The Python version used when 
developing this pipeline was Python 2.4.3 (#1, Jun 11 2009, 14:09:37). It 
should be able to run on versions newer than this. It has been tested with
Python version 2.6.6 and should work with Python version up to 2.7.x, but since
this has not been tested there is no guarantee that it will work.
 
The pipeline is distributed as a single gzipped tar-archive and is 'installed' 
simply by extracting the contents of this to a directory of your choice. It 
could be convenient to add this directory to your path, or create a symbolic
link to the program in your ~/bin directory. 
The following files are included in the archive: 
 qnrpipeline.py, 
 fluff.py, 
 readme_qnrpipeline.pdf, 
 readme_qnrpipeline.txt.
 
Here is an example of how to 'install' the pipeline:
$ tar -xf qnrpipeline-0.7030.tar.gz
$ ln -s qnrpipeline-0.7030/qnrpipeline.py ~/bin/qnrpipeline.py
In some cases it might be necessary to make the program executable by running
$ chmod +x qnrpipeline-0.7030/qnrpipeline.py


2.1 Required packages
Make sure that the following programs are installed and available on PATH:
  HMMER version 3.0 or above
  BLASTclust version 2.2.23 or above (from the BLAST package)
  MAFFT version 6.811b or above
  cdbfasta version 0.99 or above 
Please make sure that these are available and working as they should within 
your environment. To prepare the database for use with the pipeline the EMBOSS
suite toolset can be used (e.g. translation to amino-acid sequence).

Make sure that the legacy BLAST suite is installed properly, especially take
care to ensure that there is a BLASTMAT environment variable defined as this
is needed for BLASTclust to run properly. This can be tested by running
$ echo $BLASTMAT
the output you get should be a path to where the BLAST substituions matrices 
are stored. If it is not working, add the following to your ~/.bash_profile 
(make sure to enter the correct path for you system):
BLASTMAT=/path/to/blast-2.2.25/data/
export BLASTMAT



3. The pipeline


3.1 Pipeline function
The application file is 'qnrsearch.py' and can be called at the prompt if 
installed correctly and either its directory was added to the PATH-variable or
a symbolic link to the program was placed in ~/bin. Calling it without any 
arguments will produce a humble error message along with usage instructions. 
It tries to conform to standard Linux command line tool usage with optional 
flags to modify its behavior. In order to fully understand what these options 
do the way the pipeline works will now be described.


3.1.1 Pipeline layout
The pipeline can be divided into four main sub-components:
1. HMMER: hmmsearch of database using a Hidden Markov Model for the gene(s) of
   interest.
2. Identification of interesting hits through the use of a classification
   function that classifies hits based on their length of alignment to the 
   hidden Markov model together with their domain alignment bitscore according
   to some (user specifiable) function. This step also extracts the sequences,
   either from the alignments produced by hmmsearch or retrieval of sub-/
   sequences from the databases used in the search.
3. Cluster the identified sequences using BLASTclust
4. Align each cluster to study sequence similarity visually. The clusters can 
   also be aligned against reference gene sequences.
 
At each of these major steps the pipeline outputs several files so that it is 
possible to resume long runs that for some reason went wrong. It also grants 
the user the choice of running the first step once (most time consuming) and 
analysing the results by retrieving sequences of different domain scores and 
look at the clustering separately without having to redo the computationally 
intensive hmmsearch (depending on the command line switches used, this might 
not be the only very time consuming step).

In the first step there is not much settings to be made by the user, it runs
hmmsearch on the databases the user specifies and stores the output. The
second step might be more confusing since here there are a lot of things going
on "under the hood". The extraction of interesting sequences focuses on the
output from hmmsearch. The pipeline goes through the hits that hmmer found and
makes decisions on whether a hmmsearch hit is really interesting through the
use of a classification function. Here the user can decide whether to extract
the hit as it is aligned to the model by HMMER, to extract "a little more"
than what aligned in HMMER by extending the region of sequence to extract from
the source database, or finally just extract the entire source sequences. The
classification function is specially tailored to identify novel qnr-like genes
but can be set to anything the user wants, essentially. On its standard 
settings the classification function works like this:

Each fragment that matches anything in HMMER has a length, L, and a bitscore
as a measure of how well it fits the model. The fragment is classified as an
interesting sequence if the fragment bitscore is above a threshold created by
the classification function.  

It has the following parameters: 
  k - slope
  m - intercept
  d - long sequence definition
  c - threshold for long sequences

The classification function looks like this:
  classify(fragment):
   true  if fragment_score > fragment_length * k + m   if fragment_length <  d
   true  if fragment_score > c                         if fragment_length >= d
   false otherwise

So if an example fragment produces a score higher than the minimum allowed for
the current fragments length it is classified as an interesting hit and is
extracted. Such fragments are marked with * in the diagram below. Fragments
that are not classified as interesting are marked with an 'o' in the diagram
below.

       CLASSIFICATION FUNCTION
bitscore
   ^                            
   |                            
   |                            
   |  *            *            
   |          /-----------------
   |     *  /    o           o  
   |      /           o              slope: k
   |    /     o                  intercept: m
   |  /                         
   |--------------------------->
         fragment length (L)

The third step then takes the extracted sequences and clusters them using
blastclust. Here the user can modify the clustering parameters to ensure a
meaningful clustering. The fourth and last step then performs multiple
alignments of the individual clusters using mafft.


3.1.2 Command line options
The available command line options are listed below with a small description:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -n N, --numcpu=N      the number of CPUs, N, to use for programs with that
                        ability [default: N=4]
  -M MODEL, --hmm=MODEL
                        the MODEL to use in hmmsearch [default: model.hmm]
  -H, --hmmsearchonly   boolean - only run hmmsearch on the databases with the
                        model specified and save results to file. Can be
                        combined with -E. [default: False]
  -E, --extracthitsonly
                        boolean - only extract sequences from hmmsearch output
                        that classify as potential hits and then quit, writing
                        results to disk for inspection. Can be combined with
                        -H or -L. [default: False]
  -L, --clusteringonly  boolean - only perform cLustering of previously stored
                        results, requires retrieved_sequences.fasta and
                        pickled.hsseq. Can be combined with option -E.
                        [default: False]
  -A, --alignmentonly   boolean - only perform alignment of sequences within
                        previously created clusters. [default: False]
  -p PI, --percentidentity=PI
                        the percent identity, PI, with which blastclust will
                        perform the clustering, range 3-100 [default: PI=90]
  -C CT, --coveragethreshold=CT
                        the length coverage threshold, CT, used in blastclust,
                        range 0.1-0.99 [default: CT=0.25]
  -k k, --classifyK=k   modify the classification function parameter slope
                        [default: k=0.57]
  -m m, --classifyM=m   modify the classification function parameter intercept
                        [default: m=3]
  -c c, --classifyC=c   modify the classification function parameter long
                        sequence fixed cutoff [default: c=75]
  -d d, --classifyD=d   modify the classification function parameter
                        definition of long sequence (i.e. after what fragment
                        length to use fixed cutoff) [default: d=85]
  -a PATH, --addrefseq=PATH
                        add reference database with sequences in FASTA format
                        at PATH. These sequences are added before clustering
                        to improve/simplify cluster analysis. [default: not
                        used]

  Advanced options:
    Caution: some of the options might negatively affect pipeline
    performance!

    -s, --noheuristics  boolean - turn hmmsearch HEURISTICS OFF (max
                        Sensitivity) (not recommended) [default: False]
    -D, --retrdb        boolean - retrieve sequences from original database
                        and not from hmmsearch output file. Slower and not
                        'safe' for clustering or alignment. Requires an
                        cdbfasta index file for each database [default: False]
    --extendleft=EXTENDLEFT
                        An integer number of residues/nucleotides to extend
                        the HMMER aligned domain with to the left. Requires a
                        cdbfasta index file for each database [default: 0]
    --extendright=EXTENDRIGHT
                        An integer number of residues/nucleotides to extend
                        the HMMER aligned domain with to the right. Requires a
                        cdbfasta index file for each database [default: 0]

  Developer options, use at own risk!:
    Caution: some of the options might negatively affect program execution
    and/or are generally not properly tested!

    -o OUTDIR, --hmmsearch_outdir=OUTDIR
                        modify OUTput DIRectory for hmmsearch [default:
                        ./hmmsearchresults/]
    -r RESDIR, --resdir=RESDIR
                        modify the DIRectory where RESulting cluster files are
                        written [default: ./pipeline_cluster_results/]
    -R ALIGNSEQPATH, --alignseqpath=ALIGNSEQPATH
                        path to fasta file with sequences to align the
                        clusters against. It is recommended that there are not
                        more than a couple of sequences in this file. The
                        subroutine that performs the parsing is hardcoded to
                        take the five PMQR-variants! [default: False]

 
3.1.3 Pipeline output folder structure
It is recommended to perform a pipeline-run in a separate, empty, folder so 
that resulting files and folders can be created without risk of overwriting
something. The following folder structure is created upon running the entire 
pipeline, assuming the program was exectued in an empty directory called 
'pipeline':
 
pipeline/
pipeline/hmmsearchresults/
pipeline/pipeline_cluster_results/
pipeline/tmp_qnrpipe/
 
The three created folders contain, essentially, the information obtained in each 
of the three major steps of the pipeline (the fourth/alignment step, uses the 
same output folder as the third/clustering step). 
The folder 'hmmsearchresults' will contain the results from the hmmsearch of the 
databases, labeled with the run date in the filename. 
The folder 'tmp_qnrpipe' will contain the temporary files used in the 
intermediary step of the pipeline; files to enable resuming/restarting the 
clustering and multiple alignment step as well as a database file, 
'retrieved_sequences.fasta', that will contain the sequences that classified as 
potential hits according to the classification criterion. 
The last folder, 'pipeline_cluster_results', will contain two types of files 
with the results of the entire pipeline run. The files 'identified_clusters' 
and 'identified_clusters.scores' contain a flatfile representation of the 
clusters, with one cluster per line containing sequence identifiers separated by
spaces, and very similarily in the other file but here instead the sequences 
have gotten their domain score appended to the sequence identifier for simple 
identification of HMMER scores of the sequences in different clusters (note that
the reference sequences supplied with the -a switch are not present in this 
file). Also in this last folder will be several files containing the sequences 
for each cluster together with multiple alignments of each cluster.
Additionally in this folder will be clusters aligned against reference genes
specified using the command line switch -R.
 

3.2 Database preparation
The databases to be searched must be in FASTA format. To prepare the databases 
it can be convenient to use the EMBOSS toolset to convert the sequences of your
choice to FASTA format and translate all the sequences into amino-acid sequence 
using the EMBOSS tool 'transeq'. Below follows a simple example of preparing 
the database for search with the pipeline. It is assumed that you have 
downloaded or otherwise acquired a database in a suitable flat-file format (the
database, 'nucleotide', nt, from GenBank can be downloaded in FASTA format from
the NCBI FTP-server).
 
Assume the directory 'databases' contains a subdirectory 'genbank' that contains
a FASTA format file, 'nt.fasta', with the database of interest. Use for example 
EMBOSS 'transeq' to translate the database into all six frames using the 
bacterial table (number 11). The following command could be run (no
linebreaks):
[~/databases/genbank/]$ transeq -sequence nt.fasta -outseq nt.pfasta -frame=6 -table=11
 
If the user ever wishes to retrieve full-length sequences from this translated 
database (switch -D) or extend the HMMER domain alignments, an index needs to 
be created. For this the pipeline uses the program cdbfasta. Assuming that it
is installed properly, creating an index of the database file is as easy as 
running the following command:
[~/databases/genbank/]$ cdbfasta nt.pfasta
 
After successful completion of the two above commands the following three files
should be available in the ~/databases/ directory: 
nt.fasta, 
nt.pfasta,
nt.pfasta.cidx.
 

3.3 Usage example
Below follows some examples, assuming the pipeline script is made executable
and is available on the current user's PATH-variable.
 
The simplest pipeline call:
[~/pipeline/]$ qnrpipeline.py ~/databases/genbank/nt.translated.fasta
This call assumes there exists a hidden Markov model file, 'model.hmm', in 
the current directory. Take note that in real-life applications you probably
should supply at least one mandatory option; -M (to specify which model to 
use) since this makes managing what models were used when much more easy.
 
Most users might want to specify which model file should be used, and maybe 
add some reference sequences to the clustering to easier identify what 
clusters contain what kind of sequences. Then the following call could be 
used (no linebreaks):
[~/pipeline/]$ qnrpipeline.py -M ~/hmm_models/othermodel.hmm -a ~/qnrsequences/all_known_qnr_sequences.fasta ~/databases/genbank/nt.pfasta
 
For more advanced usage, assuming the pipeline has been run before, the user 
could redo the last two parts of the pipeline for analysis of a different set 
of sequences matching the Hidden Markov Model used in previous run. The 
following call will use previously stored hmmsearch result files and this time
adjust the classification parameters so that we only select sequences above 
the minimum domain score threshold of 345.0 and retrieve extended sequences 
from their original (translated) database before clustering and aligning:
[~/pipeline/]$ qnrpipeline.py -ELA -d 0 -c 345.0 --extendleft 25 --extendright 40 hmmsearchresults/* 
 
There are other options available to the user in order to specify different 
settings to BLASTclust, but they will not be discussed in detail here and the 
user is instead asked to consult the documentation provided with BLASTclust 
for detailed descriptions of their function.
 
Note that for both HMMER and BLASTclust has the ability to use standard POSIX 
threads for multiprocessor speed-ups and that the command line switch -n will
allow the user to specify the number of processors/threads available to these
parts of the pipeline, the default setting is to use four processors. For the 
programs' implementation of multi-threading, please refer to their respective 
documentation.
 
The pipeline will output a log to standard output during run time so that the
user can inspect what the running program is up to. A copy of this output is 
also written to the file 'qnrsearch.log' and stored for future reference. The
pipeline will never overwrite an old logfile but instead append the latest run
to the end of it. Please note that it WILL overwrite output files.



4. Known issues
If the user wishes to retrieve full-length sequences from the databases, 
BLASTclust might have a hard time clustering the results (i.e. it will crash)
if at least on of the sequences are complete genomes.
The pipeline currently is intended to search for Qnr-like genes and assumes 
that the five known plasmid mediated Qnr-genes are the ones to align each 
cluster against (the experimental -R option). This can be modified if the 
user desires by a simple modification in the source code in 'fluff.py',
function malign_clusters around line 1029 (and making sure that the names 
stated on this line is in the sequence identifiers in the file at 
CLUSTSEQPATH).
 
 

5. Source code
The source code is naturally available in the gzipped tar-archive since the 
Python scripts are uncompiled and the accessory file, 'fluff.py', that 
contains a host of functions important for the operation of the pipeline, is
supplied in the uncompiled state. Upon the first run of the pipeline, this 
file will compile into a byte-compiled version to improve loading times. 
Nevertheless the entire source code is available to the user for inspection 
and correction of bugs. The code should be fairly well commented and each 
function in 'fluff.py' has a docstring explaining what it does, what it 
requires as input, what output it gives and also what errors it might produce.
It can be accessed by importing it and calling their respective doc-method, 
example:
>>> import fluff
>>> print fluff.parse_blastclust.__doc__

    Parses blastclust output into a nested list structure
    
    Input:
        filename    filename of blastclust output
    Output:
        sequenceIDs list of sequence IDs with unique identifiers
                    right after the '>' symbol
    Errors:
        PathError   raised if the file does not exists
        ValueError  rasied if no unique identifiers could be found 
                    and removed
 


6. Acknowledgements
None of this would have been possible without the support of a lot of people,
helping with design choices, finding bugs or proofreading.
