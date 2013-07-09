#! /bin/bash -x

#
# Example assembly of 100bp C. elegans data set. The only argument
# this script takes is the overlap length used for the final contig assembly.
#

# We assume the data is downloaded from the SRA and converted to fastq files
# Set IN1 and IN2 to be the paths to the data on your filesystem
IN=${1##*/}
FOLDER=${1%/*}

# Parameters
SGA_BIN=sga

# Overlap parameter used for the final assembly. This is the only argument
# to the script
OL=25

# The number of threads to use
CPU=4

# To save memory, we index $D reads at a time then merge the indices together
D=4000000

# Correction k-mer value
CK=11

# The minimum k-mer coverage for the filter step. Each 27-mer
# in the reads must be seen at least this many times
COV_FILTER=1

# Overlap parameter used for FM-merge. This value must be no greater than the minimum
# overlap value you wish to try for the assembly step.
MOL=20

# Parameter for the small repeat resolution algorithm
R=10

# The number of pairs required to link two contigs into a scaffold
MIN_PAIRS=5

# The minimum length of contigs to include in a scaffold
MIN_LENGTH=200

# Distance estimate tolerance when resolving scaffold sequences
SCAFFOLD_TOLERANCE=1

# Turn off collapsing bubbles around indels
MAX_GAP_DIFF=0

cd $FOLDER
# First, preprocess the data to remove ambiguous basecalls
$SGA_BIN preprocess --pe-mode=0 -o $IN.processed --permute-ambiguous $IN

#
# Error correction
#
# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
$SGA_BIN index -a ropebwt -t $CPU --no-reverse $IN.processed


$SGA_BIN correct -a overlap -t $CPU -o $IN.corrected $IN.processed

$SGA_BIN index -a ropebwt -t $CPU $IN.corrected

# Remove exact-match duplicates and reads with low-frequency k-mers
$SGA_BIN filter -t $CPU --no-kmer-check $IN.corrected

# Merge simple, unbranched chains of vertices
$SGA_BIN fm-merge -m $MOL -t $CPU -o $IN.merged $IN.filter.pass.fa

# Build an index of the merged sequences
$SGA_BIN index -d 1000000 -t $CPU $IN.merged

# Remove any substrings that were generated from the merge process
$SGA_BIN rmdup -t $CPU -o $IN.merged.rmdup.fa $IN.merged

# Compute the structure of the string graph
$SGA_BIN overlap -m $MOL -e 0.01 -x -t $CPU $IN.merged.rmdup.fa

# Perform the contig assembly without bubble popping
$SGA_BIN assemble -m $OL -g $MAX_GAP_DIFF -r $R -o $IN.result $IN.merged.rmdup.asqg.gz

#
# Scaffolding/Paired end resolution
#
#CTGS=sga-result-contigs.fa
#GRAPH=sga-result-graph.asqg.gz

# Realign reads to the contigs
#sga-align --name aligned.pe $CTGS $IN
#echo a
# Make contig-contig distance estimates
#sga-bam2de.pl -n $MIN_PAIRS --prefix libPE aligned.pe.bam
#echo 1
# Make contig copy number estimates
#sga-astat.py -m $MIN_LENGTH aligned.pe.refsort.bam > libPE.astat
#echo 2
#$SGA_BIN scaffold -m $MIN_LENGTH --pe libPE.de -a libPE.astat -o scaffolds.n$MIN_PAIRS.scaf $CTGS
#echo z
#$SGA_BIN scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf
#echo p
