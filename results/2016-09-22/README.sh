#!/bin/bash
#
#				2016-09-23
#				----------
#
# I run here some simulations to compare the inbreeding coefficient after
# some generations of full-sib mating in Drosophila melanogaster to the
# theoretical values, determine its variance, and also to determine the mean
# and the variance of the length of inbred tracts. I use simuPOP. I simulate
# autosomes. Unfortunately, simuPOP does not support the lineage tracking
# functions in sex chromosomes.
#

if [ ! -e inbreeding_coefficient.png ]; then
   if [ ! -e noIBD_tracts.txt ]; then
      python fullsib.py > noIBD_tracts.txt
   fi
   if [ ! -e F.txt ]; then
      gawk 'BEGIN{
         print "Gen\tInd\tF"
      }(/^Gen/){
         split($1,A,/_|:/)
         F0 = length($2)
         split($2,LOCI,"")
         F = F0
         for (locus in LOCI) {
            F = F - LOCI[locus]
         }
         print A[2] "\t" A[4] "\t" F/F0
      }' noIBD_tracts.txt > F.txt
   fi
#   R --no-save < plot_F.R
fi
