#!/bin/bash
#
# This script just runs mfold on each oligonucleotide of a file that contains
# several nucleotides per row. From the output of mfold, it takes the free energy
# of the most stable fold and outputs it in the row and column corresponding to
# the oligo.

FILE=`basename $1 .txt`
FILESIZE=`wc -l $FILE | gawk '{print $1}'`
for line in `seq 1 $FILESIZE`; do
   ENERGIES=""
   for oligo in `head -n $line $FILE | tail -n 1`; do
      if echo $oligo | grep -q -P "[ACGT]*"; then
         echo ">oligo" > z$FILE.fa
         echo $oligo >> z$FILE.fa
         mfold SEQ=z$FILE.fa T=16 NA_CONC=0.05 MG_CONC=0.01 NA=DNA 1> z$FILE.log 2> z$FILE.err
         if [ -e z$FILE"_1.ct" ]; then
            ENERGIES=$ENERGIES" "`head -n 1 z$FILE"_1.ct" | gawk '{print $4}'`
         else
            ENERGIES=$ENERGIES" 0.0"
         fi
      else
         ENERGIES=$ENERGIES" NA"
      fi
      rm z$FILE*
   done
   echo $ENERGIES
done
