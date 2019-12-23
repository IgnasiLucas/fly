#!/bin/bash
#
# We have now the whole dataset available. It includes sequences from 177 libraries. I think it
# worth running dada2 outside qiime, to have more control over the process. The main reason not
# to use qiime now is that it does not let me pool the samples before denoising, which I expect
# to increase the detection power significantly.
#
# The R package dada2 is available in the qiime2-2019.10 conda environment.

DIR=('Not_used' '../../data/RUNGEN52_2019' '../../data/RUNGEN66')
SAMPLE1=(E10A E10B E10C       E11A E11B E11C       E12A E12B E12C       E14A E14B      E14D
         E15A      E15C E15D  E17A E17B E17C E17D  E19A E19B E19C E19D  E20A E20B E20C
         E22A E22B      E22D       E23B E23C E23D  E24A E24B E24C       E6A  E6B       E6D
         L10A L10B L10C L10D  L11A L11B L11C L11D  L12A L12B L12C L12D  L14A L14B L14C L14D
         L15A L15B L15C       L17A L17B L17C L17D  L19A L19B L19C L19D  L20A L20B L20CD
         L22A L22B L22C L22D  L23A           L23D  L24A L24B            L6A  L6B  L6C  L6D)

SAMPLE2=(          E22C       E23A                 E25A E25B E25C E25D  E26A E26B E26C E26D
         E27A E27B E27C       E28A E28B E28C E28D  E29A E29B E29C E29D  E30A E30B E30C E30D
         E31A E31B E31C E31D  E33A E33B E33C E33D  E35A E35B E35C E35D  E36A E36B E36C E36D
         E38A E38B E38C E38D  E39A E39B E39C E39D       L23B                      L24C
         L25A L25B L25C L25D  L26A L26B L26C L26D  L27A L27B L27C       L28A L28B L28C L28D
         L29A L29B L29C L29D  L30A L30B L30C L30D  L31A L31B L31C L31D  L33A L33B L33C L33D
         L35B L35C L35D       L36A L36B L36C L36D  L38A L38B L38C L38D  L39A L39B L39C L39D)

# After a first run of dada2, I realize that in the second sequencing batch about 50% of joined
# reads from any sample contain Ns, which are not admitted by the dada algorithm. It is unfortunate
# to remove whole reads, and so many of them, just because of probably very few Ns. I want to
# know where those Ns come from, and where they are in the reads. After looking at a couple of
# fastq files, it looks like all those Ns could be the first base of the second read, which ends
# up as the last base of the merged reads. Below, I count the number of Ns in positions of merged
# reads, starting from the last. That way, Ns in the last position are counted together, irrespectively
# of read length.

for batch in 1 2; do
   if [ ! -e Batch${batch}_N_counts.txt ]; then
      gunzip -c ${DIR[$batch]}/joined/*.extendedFrags.fastq.gz | \
      gawk -v FS="" '(NR % 4 == 2){
         for (i = 0; i < NF; i++) {
            if ($(NF - i) == "N") FREQ[i + 1]++
         }
      }END{
         for (i in FREQ) print i "\t" FREQ[i]
      }' > Batch${batch}_N_counts.txt
   fi
done

# This shows a completely different pattern of distribution of ambiguous bases between first and
# second sequencing rounds. While Ns in first sequencing-batch, joined reads are more or less
# randomly distributed along the reads, joined reads from the second batch only have Ns in last
# and third to last positions. Let's take a look at the original reverse reads:

for batch in 1 2; do
   if [ ! -e Batch${batch}_N_counts_fastq2.txt ]; then
      gunzip -c ${DIR[$batch]}/fastq/*_R2.fastq.gz | \
      gawk -v FS="" '(NR % 4 == 2){
         for (i = 1; i <= NF; i++) {
            if ($(i) == "N") FREQ[i]++
         }
      }END{
         for (i in FREQ) print i "\t" FREQ[i]
      }' > Batch${batch}_N_counts_fastq2.txt
   fi
done

# Again, the pattern is very different, suggesting different error patterns, and different
# cleaning biases between batches.
#
# The denoising analysis takes about 20 hours in Thermomix.

if [ ! -e denoising.html ]; then
   R --save -q -e "rmarkdown::render('denoising.Rmd', output_file='denoising.html')"
fi
