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

if [ ! -e F.png ] || [ ! -e NumberTracts.png ]; then
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

if [ ! -e equalRates/F.png ] || [ ! -e equalRates/NumberTracts.png ]; then
   if [ ! -e equalRates/noIBD_tracts.txt ]; then
      python fullsib_equal_rates.py > equalRates/noIBD_tracts.txt
   fi
   if [ ! -e equalRates/F.txt ]; then
      gawk -f summarize.awk equalRates/noIBD_tracts.txt > equalRates/F.txt
   fi
   R --no-save < plot_F_equal_rates.R
fi

# Now, I want to compare the standard deviation between the two recombination
# rates.

if [ ! -e SD.png ]; then
   R --no-save < plot_SD.R
fi

# According to Franklin (1977), some classic authors assumed that "after a large
# number of generations of inbreeding, the number of heterogenic segments will be
# approximately distributed as a Poisson variable, and that the distribution of
# the length of each segment (x) will tend to a negative exponential (1/a)e^(-x/a)."
# Below, I gather the lengths of the fragments found in one individual per replicate,
# in all replicates at each generation. Then, I do Q-Q plots to check the goodness
# of fit to an exponential distribution, and I estimate the parameters. Finally, I
# plot the average heterogenic fragment length with respect to the number of generations
# of full-sib mating.

if [ ! -e FragmentLength.png ] || [ ! -e QQ.png ]; then
   for i in `seq 1 19`; do
      if [ ! -e Lengths$i.txt ]; then
         gawk -v GEN=$i '($1 == GEN){
            split($10,A,",")
            for (a in A) {
                  if (A[a] != "-") print A[a]
            }
         }' F.txt | sort -nr > Lengths$i.txt
      fi
   done
   R --no-save < plot_fragments.R
fi





# References
#
# Franklin, I.R. 1977. The distribution of the proportion of the genome which is
# homozygous by descent in inbred individuals. Theoretical Population Biology 11(1):
# 60-80.
#
