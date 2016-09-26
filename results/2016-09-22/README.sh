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

if [ ! -e F.png ]; then
   if [ ! -e noIBD_tracts.txt ]; then
      python fullsib.py > noIBD_tracts.txt
   fi
   if [ ! -e F.txt ]; then
      gawk -f summarize.awk noIBD_tracts.txt > F.txt
   fi
   R --no-save < plot_F.R
fi

# With equal recombination rates in males and females:

if [ ! -d equalRates ]; then mkdir equalRates; fi

if [ ! -e equalRates/F.png ]; then
   if [ ! -e equalRates/noIBD_tracts.txt ]; then
      python fullsib_equal_rates.py > equalRates/noIBD_tracts.txt
   fi
   if [ ! -e equalRates/F.txt ]; then
      gawk -f summarize.awk equalRates/noIBD_tracts.txt > equalRates/F.txt
   fi
   R --no-save < plot_F_equal_rates.R
fi

