#!/bin/bash
#
#				2016-04-12
#				----------
#
# The question is if we can use ddRAD sequencing to massively survey the
# genetic variation in the X chromosome of Drosophila melanogaster, without
# wasting too much sequencing in the autosomes. This is a promissing strategy
# because there are chromosome X-specific sequences, and because the X
# chromosome is one of the three largest chromosomes.
#
# The strategy is to digest in silico the whole genome, and identify the
# restriction enzymes that produce the highest ratio of X to autosome DNA
# fragments.
#
# I will update a script I wrote a while ago (Cx001.plx).

DATADIR=`pwd | sed 's/results/data/'`
MAXCPUS=16
# To summarize the proportion of X fragments, use fragments of at least this
# length:
MINLENGTH=300

if [ ! -d $DATADIR ]; then
   mkdir $DATADIR
fi

# The file 'restriction' is a list of sequences recognized by restriction
# enzymes that I manually selected from REBASE. These are the sites that
# produced 3 or 5-prime overhangs, with a defined overhang sequence (with
# some exceptions) and with at least one provider. I found 94 different
# patterns. There are 8 chromosomes (or arms), including the very small
# 4 and mitochondrial. It would be more efficient to parallelize by restriction
# pattern. I will, therefore, split the patterns in several files, and
# run the files sequentially, using the maximum number of CPU cores in each.

if [ ! -e $DATADIR/restriction ]; then
   echo "The list of restriction sequences is not there."
   exit
fi

# I want to turn the list of sites in a list of perl patters suitable for
# script Cx001.plx.

if [ ! -e patterns ]; then
   date
   echo "Creating patterns."
   echo
   gawk 'BEGIN{
      X["A"] = "T"; X["C"] = "G"; X["G"] = "C"; X["T"] = "A"
      X["R"] = "Y"; X["Y"] = "R"; X["M"] = "K"; X["K"] = "M"
      X["S"] = "S"; X["W"] = "W"; X["B"] = "V"; X["D"] = "H"
      X["H"] = "D"; X["V"] = "B"; X["N"] = "N"
   }{
      REV = "";
      for (i = length($1); i >= 1; i--) {
         REV = REV X[substr($1, i, 1)]
      }
      PAT = $1
      gsub(/R/, "[GA]", PAT); gsub(/Y/, "[CT]", PAT)
      gsub(/M/, "[AC]", PAT); gsub(/K/, "[GT]", PAT)
      gsub(/S/, "[GC]", PAT); gsub(/W/, "[AT]", PAT)
      gsub(/B/, "[CGT]", PAT); gsub(/D/, "[AGT]", PAT)
      gsub(/H/, "[ACT]", PAT); gsub(/V/, "[ACG]", PAT)
      gsub(/N/, "[ACGTN]", PAT)
      if ($1 != REV) {
         PAT2 = REV
         gsub(/R/, "[GA]", PAT2); gsub(/Y/, "[CT]", PAT2)
         gsub(/M/, "[AC]", PAT2); gsub(/K/, "[GT]", PAT2)
         gsub(/S/, "[GC]", PAT2); gsub(/W/, "[AT]", PAT2)
         gsub(/B/, "[CGT]", PAT2); gsub(/D/, "[AGT]", PAT2)
         gsub(/H/, "[ACT]", PAT2); gsub(/V/, "[ACG]", PAT2)
         gsub(/N/, "[ACGTN]", PAT2)
         PAT = PAT "|" PAT2
      }
      print $1 "\t" PAT
   }' $DATADIR/restriction > patterns
fi

if [ ! -e $DATADIR/genome.fa ]; then
   for i in 2L 2R 3L 3R 4 dmel_mitochondrion_genome X Y; do
      wget ftp://ftp.ensembl.org/pub/release-84/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.chromosome.$i.fa.gz
      gunzip Drosophila_melanogaster.BDGP6.dna_sm.chromosome.$i.fa.gz
      cat Drosophila_melanogaster.BDGP6.dna_sm.chromosome.$i.fa >> $DATADIR/genome.fa
      rm  Drosophila_melanogaster.BDGP6.dna_sm.chromosome.$i.fa
   done
fi

if [ ! -e patterns1.txt ]; then
   split --suffix-length=1 --additional-suffix=".txt" --numeric-suffixes=1 --lines=$MAXCPUS patterns patterns
fi

touch checkpoints
for i in `ls -1 patterns*.txt`; do
   if ! grep -q $i checkpoints; then
      echo $i >> checkpoints
      ../../bin/Cx001.plx -g $DATADIR/genome.fa -r $i -m $MAXCPUS
   fi
done

if [ ! -e summary_$MINLENGTH.txt ]; then
   echo -e "#Numbers of fragments larger than $MINLENGTH bp" > summary_$MINLENGTH.txt
   echo -e "#Pattern\t2L\t2R\t3L\t3R\t4\tmt\tX\tY\tX/total" >> summary_$MINLENGTH.txt
   for i in `cut -f 1 patterns`; do
      gawk -v MINLEN=$MINLENGTH '($1 >= MINLEN){F[$5]++}END{for (f in F) print f "\t" F[f]}' $i/lengths | \
      gawk -v PATTERN=$i '{F[$1] = $2}END{print PATTERN "\t" F["2L"] "\t" F["2R"] "\t" F["3L"] "\t" F["3R"] "\t" F["4"] "\t" F["dmel_mitochondrion_genome"] "\t" F["X"] "\t" F["Y"] "\t" F["X"] / (F["2L"] + F["2R"] + F["3L"] + F["3R"] + F["4"] + F["dmel_mitochondrion_genome"] + F["X"] + F["Y"])}' >> summary_$MINLENGTH.txt
   done
fi

# Conclusions
# -----------
#
# Using a single enzyme, instead of two, we get an idea of the maximum number of fragments
# that we can produce per chromosome. While a combination of two enzymes could increase the
# proportion of X-specific fragments, it would reduce the total number of fragments.
#
# The number of fragments larger than 500 bp that can be obtained from the X chromosome range between 152
# (pattern GGCGCGCC) and 16720 (pattern ACGT). The average proportion of X fragments is
# 17.4%, and its range is 15.0 to 27.2%.
#
# These numbers are not as promissing as I hoped. I conclude that we would have to sequence
# a lot, and waste a lot of sequences, to survey the variation in the X chromosome with enough
# detail.
#
# I am still wanting to try a GBS experiment. I would aim at no more than 3000 fragments in
# chromosome X. I require the overhang not to include the whole restriction site, so that I can
# design adapters that do not reproduce the restriction site, and use afterwards the same enzyme
# to get rid of chimeras. One interesting candidate is NspI, which recognizes RCATGY, leaves a
# CATG-3' overhang, and produces 19043 fragments between 250 and 500 bp, including 3352 from the
# X chromosome. Another reasonable candidate is HaeII, cutting on RGCGCY, leaving a GCGC3' overhang
# and producing 11455 genomic fragments between 250 and 500, including 2339 from chromosome X.
#
# In order to try another enzyme, I want to know how many fragments of each type are produced
# by some frequent cutters with degenerate recognition sites. In particular, I am interested in
# pattern GDGCHC, recognized by Bsp1286I, and leaving DGCH-3' overhang
#
# GDGCHC	G[AGT]GC[ACT]C
#
# Such a degenerate recognition site corresponds to 6 different sequences (see below), only three
# of which are palindromic. I wondered whether all 6 sequences would be cut by the enzyme or if it
# would only cut the palindromic ones. I digested in silico the genome of the bacteriophage M13mp7
# and reproduced the pattern of 4 fragments experimentally obtained with the isoschizomer BmyI,
# as reported by  Wagner et al. (1990; Nucl. Acids Res. 18[10]:3088). The 4 cut sites in M13mp7 are
# all different, and only one is palindromic (not shown).
#
# How many types of fragments would Bsp1286I generate? The 3 non-palindromic sequences can be in
# 2 different orientations, relative to the site in the opposite extreme. Thus, they count as double
# number of sequences. Regardless of the order in which they appear, and allowing the same sequence
# in both ends, this amounts to 9 + (9 over 2) = 9 + 9 * 8 / 2 = 45 fragment types.

#if [ ! -e patterns ]; then
#   echo -e "GDGCHC1\tGAGCAC|GTGCTC"  > patterns2 # AGCA-3'/TGCT-3' overhangs.
#   echo -e "GDGCHC2\tGAGCCC|GGGCTC" >> patterns2
#   echo -e "GDGCHC3\tGAGCTC"        >> patterns2
#   echo -e "GDGCHC4\tGGGCAC|GTGCCC" >> patterns2
#   echo -e "GDGCHC5\tGGGCCC"        >> patterns2
##  echo -e "GDGCHC6\tGGGCTC|GAGCCC" >> patterns2 # same as GDGCHC2
#   echo -e "GDGCHC7\tGTGCAC"        >> patterns2
##  echo -e "GDGCHC8\tGTGCCC|GGGCAC" >> patterns2 # same as GDGCHC4
##  echo -e "GDGCHC9\tGTGCTC|GAGCAC" >> patterns2 # same as GDGCHC1
#fi

