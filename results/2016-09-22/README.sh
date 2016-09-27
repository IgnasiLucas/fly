#!/bin/bash
#
#				2016-09-23
#				----------
#
# I run here some simulations to compare the inbreeding coefficient after
# some generations of full-sib mating in Drosophila melanogaster to the
# theoretical values, determine its variance, and also to determine the mean
# and the variance of the length of inbred tracts. I use simuPOP. I simulate
# autosomes. SimuPOP does not support the lineage tracking functions in sex
# chromosomes.
#
# Unfortunately, the parameters of the simulations are hard coded in the python
# script. Among them, the number of markers along the chromosome must be set
# together with the distance between markers. The approach here is to monitor
# markers at equal genetic distances. A fly autosome is around 100 cM, meaning
# that one cross over is expected per meiosis, I think. In simuPOP the recombination
# rate between markers is the probability of a recombination happening there.
# I set the rate to the reciprocal of the number of markers. I understand that
# a recombination rate of 0.01 is equivalent to 1 cM. The average recombination
# intensity in Drosophila is around 2 cM/Mb, in autosomes. Thus, on average a
# 10 cM fragment is around 5 Mb.
#


if [ ! -e F.png ] || [ ! -e NumberTracts.png ]; then
   if [ ! -e noIBD_tracts.txt ]; then
      python fullsib.py > noIBD_tracts.txt
   fi
   if [ ! -d equalRates ]; then mkdir equalRates; fi
   if [ ! -e equalRates/noIBD_tracts.txt ]; then
      python fullsib_equal_rates.py > equalRates/noIBD_tracts.txt
   fi
   if [ ! -e F.txt ]; then
      gawk -f summarize.awk noIBD_tracts.txt > F.txt
   fi
   if [ ! -e equalRates/F.txt ]; then
      gawk -f summarize.awk equalRates/noIBD_tracts.txt > equalRates/F.txt
   fi
   R --no-save < plot_F.R
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

if [ ! -e FragmentLength.png ] || [ ! -e QQ.png ] || [ ! -e equalRates/QQ.png ]; then
   for i in `seq 1 19`; do
      if [ ! -e Lengths$i.txt ]; then
         gawk -v GEN=$i '($1 == GEN){
            split($10,A,",")
            for (a in A) {
               if (A[a] != "-") print A[a]
            }
         }' F.txt | sort -nr > Lengths$i.txt
      fi
      if [ ! -e equalRates/Lengths$i.txt ]; then
         gawk -v GEN=$i '($1 == GEN){
            split($10,A,",")
            for (a in A) {
               if (A[a] != "-") print A[a]
            }
         }' equalRates/F.txt > equalRates/Lengths$i.txt
      fi
   done
   R --no-save < plot_fragments.R
fi

# CONCLUSIONS
# -----------
#
# Figure F.png gives me confidence that the simulations are well
# done, since the average inbreeding coefficient is almost identical to the theoretical
# expectation. Note that whether males undergo recombination or not does not affect the average
# inbreeding coefficient. Indeed, recombination rate is not expected to affect the
# average inbreeding coefficient, but its variance (see Franklin, 1977). Comparing the
# two plots, the difference in the variance of F between absence and presence of male
# recombination is noticeable. Figure SD.png shows that the standard deviation of F along
# the generations of inbreeding is higher in the absence of recombination in males (red
# line) than in its presence.
#
# Unfortunately, simuPOP does not support the tracking of allele lineages in sex chromosomes.
# However, it is quite clear what we can expect. Because the X chromosome spends more time
# in females, it must experience more cross overs than autosomes (in D. melanogaster). Part
# of the excess variance in F in the absence of male recombination must be due to the variance
# of the number of generations a chromosome spent in females. In only a handful of generations
# before inbreeding is complete, chromosomes must differ substantially in the number of
# generations they spent in either males or females. For the X chromosome, however, this variance
# is lower, because it cannot spend two consecutive generations in males. Thus, I expect the
# X chromosome along a full-sib mating programme in D. melanogaster to accumulate IBD
# tracts in a less variable pace than autosomes.
#
# In figures NumberTracts.png and equalRates/NumberTracts.png, the colours indicate the number
# of segments of a chromosome that are not identical by descent yet, along the generations
# of inbreeding. The clearest yellow corresponds to 0 fragments; that is, completely homozygous
# chromosomes. It is striking to see how much higher the proportion of fully IBD chromosomes
# is in the absence than in the presence of recombination in males. An implication of this is
# that in the absence of male recombination, autosomes run a much higher risk of becoming
# completely homozygous in just a few generations of inbreeding (20% chance by the 6th generation
# of full-sib mating) than in the presence of male recombination (5% chance).
#
# It would be interesting to know what the situation is for the X chromosome. Since it must
# spend not 50%, but 66.7% of its time in females, the variance in its inbreeding coefficient must
# be even lower than that observed in autosomes under equal recombination rates between males
# and females. Thus, along the inbreeding experiment, I expect some lines to turn accidentally
# homozygous for one or even both autosomes before losing the heterozygosity in the X chromosome.
# However, the presence of recessive lethals will make these cases look less frequent among
# the surviving lines.
#
#
# References
#
# Franklin, I.R. 1977. The distribution of the proportion of the genome which is
# homozygous by descent in inbred individuals. Theoretical Population Biology 11(1):
# 60-80.
#
