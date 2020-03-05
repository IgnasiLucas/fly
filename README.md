# Fly

This is a collaboration with Pau Carazo and Zahida Sultanova. We are comparing the gut
microbiome of several Drosophila melanogaster isolines, at two time points in their
life cycle. We sequenced 16S rRNA gene amplicons and we are analyzing them with dada2.

------------------------------------------------------------------------------------

## 2020-02-11
I use dada2 to assign taxonomy to the amplicon sequence variants. Among 2892 sequences,
only 16 are annotated to the species level. But 2764 (96%) are annotated to the genus,
which is quite good. There are 24 eukaryotic sequences, the annotation of which failed,
probably because of the database not being appropriate. The rest are bacteria. This is
The distribution of bacterial phyla (not considering abundances, but just number of variants):

- Proteobacteria:        2303
- Firmicutes:             382
- Actinobacteria:          92
- Bacteroidetes:           44
- Cyanobacteria:           30
- Verrucomicrobia:          4
- Deinococcus-Thermus:      3
- Fusobacteria:             3
- Acidobacteria:            2
- Planctomycetes:           2
- Chloroflexi:              1
- Gemmatimonadetes:         1
- Patescibacteria:          1

See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-02-11/taxonomy.html).

## 2019-12-12
I run dada2 on the whole dataset, pooling samples before calling the unique variants.
It took around 20 hours. Unfortunately, I had to trim reads and remove those that still
included Ns. The second sequencing batch produced a much more variable depth among
samples; isolines 27 and 35 are under-sequenced early in life, while isoline 29 is
badly over-sequenced late in life. See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-12-12/denoising.html)
and the alignment of the sequences [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-12-12/alignment.html).

## 2019-12-11
The second batch of sequences arrived. Three of the 28 combinations of isoline and time
point lack replicates. There are large differences in number of reads and in merging
success rate among samples. The dominant merged lengths are the same as in the first
batch. See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-12-11/explore.html).
But, again, be aware that a table is not rendering properly in that preview.

## 2019-11-25
I reproduce the plots and tables offered by qiime's visualization of the table of
'features' (a.k.a. exact sequence variants). I also produced a heatmap to visualize
the whole abundance table. See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-11-25/ExploreData.html)

Unfortunately, the tables from the *reactable* package do not work in the link above.
Please, download the 2019-11-25/ExploreData.html file and open it in your browser.

There are 2472 exact sequence variants. The most frequent abundance is also the lowest,
2 observations in one sample. I suspect that a de-noising procedure with pooled samples
would detect many more, low-frequency variants.

## 2019-11-14
I imported the sequence data into qiime and de-noised the reads. I did it starting from
both, the already merged reads, and the clean but not merged yet read pairs. This last
way of processing the data is not necessary, but I was interested in checking how it
works, and to evaluate the merged dataset offered by the sequencing services. I learned
that dada2 trimming reads before merging and discarding reads shorter than the truncating
point is bad practice, based on a flawed assumption: that fragments from
which reads originated are all of very similar lengths. This contrasts with the merged
data available, with two high-frequency lengths (440, 465). In this situation, it seems
that dada2 default settings in qiime2 can introduce a bias in composition. See the
report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-11-14/merging_optimization.html)

## 2019-11-11
We got the fastq files from the first batch of amplicon sequencing. The provider run some
preliminary analysis, and delivered a report on the results of the quality control process
for every sample. Here, I pool all samples together to look at the numbers of reads and
their length distributions along the pipeline. You can see a report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-11-11/explore.html)

## 2019-06-11
k-mer analysis of some preliminar fly gut microbiota metagenomics data. My goal is to
assess data quality and levels of diversity before assembly. You can see rendered
versions of the reports in the following links:

* [Analysis of sequence redundancy](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-06-11/Sequence_Redundancy.html)
* [GC and k-mer frequency heatmaps](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-06-11/GCP_heatmaps.html)
* [k-mer frequency-based comparison among samples](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2019-06-11/Comparisons.html)

## 2016-12-15
Use freebayes with the four fastq files that contain unambiguously identified reads,
in order to identify diagnostic SNPs. But I interrumpted this analysis here, because
the real experiment of which this was a pilot got cancelled, for strategical reasons.

## 2016-12-09
I want to call variants and estimate the potential heterozygosity among lines. I want to
use both mpileup and freebayes, and maybe GATK as well.

## 2016-12-08
I run step 1 of pyrad (demultiplexing), to compare the results with those of sabre, and
to get an idea of how many more reads we can save from the undetermined fastq files.
I also generate the histograms of merged read lengths.

## 2016-12-03
I use bedtools to compare the bam files among the 4 samples. It is clear that
sample 1, processed with only 12 PCR cycles, has a better distribution of coverage
per site. All samples have more than 12000 sites covered at least 6 times, and
more than 9000 of them are covered at least 6 times in all 4 samples.

## 2016-11-28
Merged paired reads, trimmed adapters, and mapped them to the reference genome.


## 2016-11-25
Demultiplexed the fastq files provided by the sequencing center. I am keeping track of the
contradictions between the indices and the in-line barcode.


## 2016-10-04
The goal here is to calculate the composition of the ligation reactions, in order to have
a 10 fold excess of adapters to fragment ends

## 2016-09-22
I use simuPop to simulate 20 generations of full-sib mating, and estimate the effect of
inbreeding on the tracts of an autosome that are not identical by descent yet. The simulations
reproduce the expected increase of inbreeding coeficient (F). However, the variance of F
is high in general, and even higher in the absence of recombination in the male germline of
Drosophila melanogaster. As a consequence of this variance, we expect almost 50% of flies
inbred by full-sib mating for only 10 generations to be already completely homozygous.

## 2016-06-09
Design of 12 adapters with in-line 8-nucloetide codewords, for use with restriction enzyme
NspI.

## 2016-04-12

The original goal was to check if the presence of X-specific sequences in Drosophila,
involved in the gene-dose compensation system, allowed for an enrichment of X-chromosome
fragments using restriction enzymes. The conclusion was negative. But I decided to keep
pursuing a GBS experiment in Drosophila, and used the results to choose the most appropriate
enzymes: NspI, and HaeII. I do not think we need more than 3000 fragments in the X chromosome

Assuming a yield of 20 million reads, and 20000 different genomic fragments, we could genotype
20 flies at 50X coverage in one MiSeq run.
