#!/bin/bash
#
# This script just runs mfold on each oligonucleotide of a file that contains
# several nucleotides per row. From the output of mfold, it takes the free energy
# of the most stable fold and outputs it in the row and column corresponding to
# the oligo.

FILE=`basename $1 .txt`
FILESIZE=`wc -l $FILE.txt | gawk '{print $1}'`
for line in `seq 1 $FILESIZE`; do
   ENERGIES=""
   for oligo in `head -n $line $FILE.txt | tail -n 1`; do
      if echo $oligo | grep -q -P "^[ACGT]*$"; then
         echo ">oligo" > z$FILE.fa
         echo $oligo >> z$FILE.fa
         #mfold writes a few lines to /dev/tty. They are caught neither by stderr nor by stdout.
         mfold SEQ=z$FILE.fa T=16 NA_CONC=0.05 MG_CONC=0.01 NA=DNA MAX=10 1> z$FILE.log 2> z$FILE.err
         if [ -e z$FILE".fa_1.ct" ]; then
            ENERGIES=$ENERGIES"\t"`head -n 1 z$FILE".fa_1.ct" | gawk '{print $4}'`
         else
            ENERGIES=$ENERGIES"\t0.0"
         fi
         rm z$FILE*
      else
         ENERGIES=$ENERGIES"\tNA"
      fi
   done
   echo -e $ENERGIES
done
