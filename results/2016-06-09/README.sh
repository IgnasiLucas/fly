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
        gawk -v COMP=1 -v OFS="\t" -f ../../bin/Cx004.awk codes.top > codes.bot
fi

# [Note: I already checked that the minimum Hamming distance is 4].
