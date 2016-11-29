Fly
===

This is a collaboration with Pau Carazo and Zahida Sultanova. We need to genotype
some flies at several loci along the X chromosome. The SNP arrays commercially available
are too expensive, and would only make sense if we had to genotype a very large number
of samples. There are other options, that I am not familiar with yet. Since I am interested
in updating the GBS protocol, I will prepare a pilot experiment.

As usual, below I explain the goal of the analyses run in each folder, in reverse
chronological order.

2016-11-28
==========
Merged paired reads, trimmed adapters, and mapped them to the reference genome.


2016-11-25
==========
Demultiplexed the fastq files provided by the sequencing center. I am keeping track of the
contradictions between the indices and the in-line barcode.


2016-10-04
==========
The goal here is to calculate the composition of the ligation reactions, in order to have
a 10 fold excess of adapters to fragment ends

2016-09-22
==========
I use simuPop to simulate 20 generations of full-sib mating, and estimate the effect of
inbreeding on the tracts of an autosome that are not identical by descent yet. The simulations
reproduce the expected increase of inbreeding coeficient (F). However, the variance of F
is high in general, and even higher in the absence of recombination in the male germline of
Drosophila melanogaster. As a consequence of this variance, we expect almost 50% of flies
inbred by full-sib mating for only 10 generations to be already completely homozygous.

2016-06-09
==========
Design of 12 adapters with in-line 8-nucloetide codewords, for use with restriction enzyme
NspI.

2016-04-12
==========

The original goal was to check if the presence of X-specific sequences in Drosophila,
involved in the gene-dose compensation system, allowed for an enrichment of X-chromosome
fragments using restriction enzymes. The conclusion was negative. But I decided to keep
pursuing a GBS experiment in Drosophila, and used the results to choose the most appropriate
enzymes: NspI, and HaeII. I do not think we need more than 3000 fragments in the X chromosome

Assuming a yield of 20 million reads, and 20000 different genomic fragments, we could genotype
20 flies at 50X coverage in one MiSeq run.
