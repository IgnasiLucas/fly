Fly
===

This is a collaboration with Pau Carazo and Zahida Sultanova. We need to genotype
some flies at several loci along the X chromosome. The SNP arrays commercially available
are too expensive, and would only make sense if we had to genotype a very large number
of samples. There are other options, that I am not familiar with yet. Since I am interested
in updating the GBS protocol, I will prepare a pilot experiment.

As usual, below I explain the goal of the analyses run in each folder, in reverse
chronological order.

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
