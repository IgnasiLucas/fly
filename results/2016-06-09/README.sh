#!/bin/bash
#
#				2016-06-09
#				----------
#
# I want to design the adapters for a GBS experiment. Because the genomes
# will be fragmented with only one frequent cutter (either HaeII or NspI)
# I need divergent, Y-shaped adapters. I have already designed and tested
# a set of 12 codewords of sizes between 5 and 8 bases. I can combine those
# with indices, added during the PCR amplification.
#

if [ ! -e illumina-adapter-sequences_1000000002694-01.pdf ]; then
   wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf
fi

#
# This is the Illumina Truseq universal adapter, with an index.
#
#         5-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3
#            | ||                      |||      |        ||||||||||||
#   3-GTTCGTCTTCTGCCGTATGCTCTAGCACTACACTGACCTCAAGTCTGCACACGAGAAGGCTAG-5
#                             ------
#                             index1
#
# This is what Andolfatto et al. (2011) used (where Ns are the barcode):
#
# P2      5-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA-3
#                                  5-ACACTCTTTCCCTACACGACGCTCTTCCGANNNNNN
#                                    |  |      |         ||||||||||||||||
#                                 3-GTTCGTCTTCTGCCGTATGCTCGAGAAGGCTnnnnnnAT-5
#
#
# And this is one of the adapters used by Peterson et al. (2012):
#
#                                                                     brcd.
#         5-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACG-3          -----
#                                  5-ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATG-3
#                                    ||||||||||||||||||||||||||||||||||||||
#                                   3TGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGACGTACTTAA5
#
#
# The 'Illumina Adapter Sequences' document, from February 2016, confirms that the flow-cell
# attachment sequences did not change. I also have confirmation from the sequencing services
# that their MiSeq machine is sensitive to the nucleotide composition in the first few cycles.
# Thus, I plan to use the sets of balanced codewords that I had prepared for Culex pipiens,
# and this time in combination with an index. The convenience of adding an index justifies the
# use of adapters shorter than the TruSeq version, and complemented later during the amplification
# PCR. An additional reason to use short adapters may be the higher efficiency of the oligo
# synthesis for shorter oligos. In addition, I will keep applying the reduction of the Illumina
# TruSeq adapter and its sequencing primer to make some room for the codeword (Andolfatto et al.
# 2011), even though it may not be necessary (Peterson et al. [2012] do not apply it), and it
# forces me to use customized sequencing primers.
#
# I ordered two enzymes, NspI and HaeII, with different cut sites. I will design both sets of
# adapters, but eventually order only one (for NspI).
#
# +========+=========+==============+===========+
# | Enzyme |   Site  | 250-500 tot. | 250-500 X |
# +--------+---------+--------------+-----------+
# | NspI   | RCATG'Y |        19043 |      3352 |
# |        | Y'GTACR |              |           |
# |        |         |              |           |
# | HaeII  | RGCGC'Y |        11455 |      2339 |
# |        | Y'CGCGR |              |           |
# +========+========================+===========+
#
# NspI-adapters
# =============
#
# P2      5-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA-3
#                                  5-ACACTCTTTCCCTACACGACGCTCTTCCGAXXXXXXXYCATG-3
#                                      |||      |        ||||||||||||||||||
#                                 3-CACTGACCTCAAGTCTGCACACGAGAAGGCTxxxxxxxR-5
#
# The adapters are composed of a top and a bottom oligo. The top oligo has 30 constant
# nucleotides, followed by an 8 nucleotides codeword, a sufix of 0 to 3 nucleotides and
# the overhang that matxes that of the digested fragments. In order to avoid the regeneration
# of the restriction site upon ligation of the adapters to the genomic fragments, the last
# nuclotide before the overhang must be C or T. This will allow the secondary restriction
# of chimeric fragments, or the simultaneous digestion and ligation.
#
# It is not difficult to generate a set of 4000 codes of size 12, length 8 and balanced
# composition that guarantee a minimum Hamming distance of 4:

if [ ! -e codes.top ] || [ ! -e codes.bot ]; then
        ../../bin/Cx003.plx | grep -v -P "AAA|CCC|GGG|TTT" > codes.top
fi

# For simplicity, I balanced each set of 4 consecutive codewords in each code.
# [Note: I already checked that the minimum Hamming distance is 4].
# The 8-nucleotide codewords having balanced composition should be enough for the MiSeq
# to calibrate correctly during those first cycles of sequencing. However, it is not
# difficult to add a linker sequence between the codeword and the overhang in order
# to stagger the overhangs and extend the balanced compositions of nucleotides at least
# over the following four sequencing cycles. This is particularly convenient for NspI
# overhangs (CATG). The linker and overhang could be the following.
#
#       CATG
#       TCATG
#       GTCATG
#       AGTCATG
#
# Each linker-overhang sequence should be added to 3 of the codewords of a code, because the
# code has 12 codewords. The question is to distribute these four possible sufixes among the
# codewords in such a way that:
#
#   1. The shortest sufix is added only to codewords ending in T or C.
#   2. The potential secondary structure of DNA is minimized among the resulting set of
#      oligos (including top and bottom).
#   3. And the complementarity among all the oligos is limited to the expected stretch
#      between top and bottom of a pair.
#
# Any distribution of the 4 sufixes above among the 12 codewords that produce 3 oligos of
# each size would, in principle, be valid. There are (12 over 3) possible ways to choose
# the 3 codewords of one size class, (9 over 3) possible ways to choose the 3 or the following
# size class, and (6 over 3) possible ways to choose the 3 codewords of the third size class,
# the remaining 3 being necessarily of the 4th size class. Thus, this would amount to evaluating
# 369600 distributions of the sufixes among codewords per code, times 4036 codes = 1491705600
# possible oligo sets.
#
# However, I can constrain the distribution to have the four possible sizes in
# each set of four consecutive codewords, so that I keep each set balanced, which allows me
# to use only 4 or 8, instead of necessarily 12 oligos at any time. With this constraint,
# there are 4! x 4! x 4! = 13824 possible size distributions per code, which should be easier to
# evaluate. They are even fewer valid ones, if we remove the size distributions that reproduce
# the restriction site. Indeed if I add the shortest suffix only to codewords ending in T
# (one per set of 4 consecutive codewords), I get only 3! x 3! x 3! = 216 combinations per code,
# or 871776 possible sets of oligonucleotides. Which one is best?

if [ ! -e oligos ]; then
   gawk 'BEGIN{
      TOPPREFIX = "ACACTCTTTCCCTACACGACGCTCTTCCGA"
      TOPSUFIX[0] = "CATG"; TOPSUFIX[1] = "TCATG"; TOPSUFIX[2] = "GTCATG"; TOPSUFIX[3] = "AGTCATG"
      BOTPREFIX[0] = ""; BOTPREFIX[1] = "A"; BOTPREFIX[2] = "AC"; BOTPREFIX[3] = "ACT"
      BOTSUFIX = "TCGGAAGAGCACACGTCTGAACTCCAGTCAC"
      COMP["A"] = "T"; COMP["C"] = "G"; COMP["G"] = "C"; COMP["T"] = "A"
   }
   function revcomp(SEQ,   REV){
      REV = ""
      for (i = length(SEQ); i > 0; i--) {
         REV = REV COMP[substr(SEQ, i, 1)]
      }
      return REV
   }
   {
      for (i = 1; i <= NF; i++) {
         CODEWORD[i] = $i
      }
      for (SET = 0; SET <= 2; SET++) {
         for (CWINSET = 1; CWINSET <= 4; CWINSET++) {
            CODEWORD[substr($(SET * 4 + CWINSET), 8, 1)] = $(SET * 4 + CWINSET)
         }
         T = 0
         i = 0
         for (A = 1; A <= 3; A++) {
            for (C = 1; C <= 3; C++) {
               if (C != A) {
                  for (G = 1; G <= 3; G++) {
                     if ((G != A) && (G != C)) {
                        i++
                        TOPOLIGO[SET, i, "A"] = TOPPREFIX CODEWORD["A"] TOPSUFIX[A]
                        TOPOLIGO[SET, i, "C"] = TOPPREFIX CODEWORD["C"] TOPSUFIX[C]
                        TOPOLIGO[SET, i, "G"] = TOPPREFIX CODEWORD["G"] TOPSUFIX[G]
                        TOPOLIGO[SET, i, "T"] = TOPPREFIX CODEWORD["T"] TOPSUFIX[T]
                        BOTOLIGO[SET, i, "A"] = BOTPREFIX[A] revcomp(CODEWORD["A"]) BOTSUFIX
                        BOTOLIGO[SET, i, "C"] = BOTPREFIX[C] revcomp(CODEWORD["C"]) BOTSUFIX
                        BOTOLIGO[SET, i, "G"] = BOTPREFIX[G] revcomp(CODEWORD["G"]) BOTSUFIX
                        BOTOLIGO[SET, i, "T"] = BOTPREFIX[T] revcomp(CODEWORD["T"]) BOTSUFIX
                     }
                  }
               }
            }
         }
      }
      n = 0
      for (i = 1; i <= 6; i++) {
         for (j = 1; j <= 6; j++) {
            for (k = 1; k <= 6; k++) {
               n++
               CODE = sprintf("CODE%04u.%03u", NR, n)
               CODE = CODE "\t" TOPOLIGO[0, i, "A"] "\t" TOPOLIGO[0, i, "C"] "\t" TOPOLIGO[0, i, "G"] "\t" TOPOLIGO[0, i, "T"]
               CODE = CODE "\t" TOPOLIGO[1, j, "A"] "\t" TOPOLIGO[1, j, "C"] "\t" TOPOLIGO[1, j, "G"] "\t" TOPOLIGO[1, j, "T"]
               CODE = CODE "\t" TOPOLIGO[2, k, "A"] "\t" TOPOLIGO[2, k, "C"] "\t" TOPOLIGO[2, k, "G"] "\t" TOPOLIGO[2, k, "T"]
               CODE = CODE "\t" BOTOLIGO[0, i, "A"] "\t" BOTOLIGO[0, i, "C"] "\t" BOTOLIGO[0, i, "G"] "\t" BOTOLIGO[0, i, "T"]
               CODE = CODE "\t" BOTOLIGO[1, j, "A"] "\t" BOTOLIGO[1, j, "C"] "\t" BOTOLIGO[1, j, "G"] "\t" BOTOLIGO[1, j, "T"]
               CODE = CODE "\t" BOTOLIGO[2, k, "A"] "\t" BOTOLIGO[2, k, "C"] "\t" BOTOLIGO[2, k, "G"] "\t" BOTOLIGO[2, k, "T"]
               print CODE >"oligos"
            }
         }
      }
   }' codes.top
fi

PROCESSES=24
rm z*
split -d --additional-suffix=.txt -n l/$PROCESSES oligos zoligos

# Add control, and call the following in a script, to parallelize
for file in `basename -s .txt zoligos*.txt`; do
   if [ ! -e $file.mfold ] && [ ! -e oligosdata ]; then
      ../../bin/runmfold.sh $file.txt > $file.mfold &
   fi
done
wait


# The files *.mfold contain the energies of the most stable secondary structures
# of all the oligos. They are usually negative numbers. The higher the number, the
# less stable the structure is, and the better for the purpose of the oligos. Both
# bottom and top oligos must be maximized.
#
# In addition, I need to check the self-alignments. I have checked that mofld does not
# consider dimers. That is a limitation, since oligos can hybridize with one another,
# in addition to hybridizing within the same molecule. The problem of dimer prediction
# is one of DNA alignment, although phrased differently. Alignments search for similarity,
# while hybridization prediction should search for complementarity. However, complementarity
# is but similarity to the reverse complementary strand. Thus, only alignments between
# the two different strands are signs of potential hybridization.
#
# I can use exonerate, with a low minimum score, in order to report non-trivial, suboptimal
# alignments. The complementarity between the top and the bottom strands of the adapters
# extends through at least 18 base pairs, which seems to give a score of 90 in exonerate,
# with default penalties. Although melting temperatures, rather than alignment scores
# should be used, I think that any local alignment with score < 20 will not compete
# strongly with the expected hybridization.

for file in `basename -s .txt zoligos*.txt`; do
   if [ ! -e $file.exonerate ] && [ ! -e oligosdata ]; then
      ../../bin/runexonerate.sh $file.txt > $file.exonerate &
   fi
done
wait

if [ ! -e oligosdata ]; then
   cat zoligos*.mfold > oligos.mfold
   cat zoligos*.exonerate > oligos.exonerate
   paste oligos.mfold oligos.exonerate | cut -f 2-25,27- > oligosdata
   rm oligos.mfold oligos.exonerate zoligos*
fi

# The file oligosdata should be a matrix of 871776 rows and 48 columns. Each row corresponds to
# a possible set of 12 top and 12 bottom oligos. The 48 columns are the 24 free energies of
# the best folds of each oligo, and the 24 scores of the best alignments to self-reverse-
# complements. Now the question is to find the row with maximum free energies and minimum alignment
# scores.

if [ ! -e trace ]; then
   gawk 'BEGIN{
      AVEENERGY = -5.0
      MINENERGY = -5.0
      AVEDIMERS	= 100
      MAXDIMERS = 100
   }{
      SUMENERGY = 0
      SUMDIMERS = 0
      THISMINENERGY = 0
      THISMAXDIMERS = 0
      for (i = 1; i <= 24; i++) {
         SUMENERGY += $i
         if ($i < THISMINENERGY) THISMINENERGY = $i
      }
      for (i = 25; i <= 48; i++) {
         SUMDIMERS += $i
         if ($i > THISMAXDIMERS) THISMAXDIMERS = $i
      }
      THISAVEENERGY = SUMENERGY / 24
      THISAVEDIMERS = SUMDIMERS / 24
      if ((THISAVEENERGY >= AVEENERGY) && (THISMINENERGY >= MINENERGY) && (THISAVEDIMERS <= AVEDIMERS) && (THISMAXDIMERS <= MAXDIMERS)) {
         AVEENERGY = THISAVEENERGY
         MINENERGY = THISMINENERGY
         AVEDIMERS = THISAVEDIMERS
         MAXDIMERS = THISMAXDIMERS
         BESTOLIGO = $1
         print BESTOLIGO "\t" AVEENERGY "\t" MINENERGY "\t" AVEDIMERS "\t" MAXDIMERS
      }
   }' oligosdata > trace
fi

# The last row of file trace should have the name of the best set of oligos.
