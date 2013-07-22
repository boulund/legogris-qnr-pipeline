#!/bin/bash
echo -e "File\tLine count\tLongest line"
echo -e "------------------------------------------------------"
for f in `ls -d */ | awk '{ print $1 }'`
do
    file="${f}fragments_passed.nfa.result-contigs.fa"
    echo -e "$f\t`wc -l $file | awk '{print $1}'`\t`wc -L $file | awk '{print $1}'`"
done
