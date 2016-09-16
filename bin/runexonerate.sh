#!/bin/bash

FILE=`basename $1 .txt`
FILESIZE=`wc -l $FILE.txt | gawk '{print $1}'`
for line in `seq 1 $FILESIZE`; do
   SCORES=""
   for oligo in `head -n $line $FILE.txt | tail -n 1`; do
      # Here we go through every oligo in a row of the original file.
      # But the first word in every line is an identifier, which will
      # not be processed with exonerate.
      if echo $oligo | grep -q -P "^[ACGT]*$"; then
         echo ">oligo" > z$FILE.fa
         echo $oligo  >> z$FILE.fa
         exonerate --query  z$FILE.fa \
                   --target z$FILE.fa \
                   --model affine:local \
                   --exhaustive \
                   --percent 0 \
                   --score 20 \
                   --showalignment false 2> /dev/null | \
         gawk -v S=19 '((/^vulgar/) && (S == 19) && ($5 != $9)){
            S = $10; print S; nextfile}' > z$FILE.score
         SCORES=$SCORES"\t"`cat z$FILE.score`
      else
         SCORES=$SCORES"\tNA"
      fi
      # Unfortunately, the string SCORES will start with a tab.
   done
   echo -e $SCORES
   rm z$FILE*
done
