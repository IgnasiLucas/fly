# Fly

This is a collaboration with Pau Carazo and Zahida Sultanova. We are comparing the gut
microbiome of several Drosophila melanogaster isolines, at two time points in their
life cycle. We sequenced 16S rRNA gene amplicons and we are analyzing them with dada2.

------------------------------------------------------------------------------------

## 2020-05-01
Amparo suggested to use the linear discriminant analysis, implemented
in the Galaxy tool Lefse, to identify biomarkers of age, or other categorical
variables. In this folder I prepare the data and I link the results in the Galaxy
session. The report is [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-05-01/lefse.html)
and the Galaxy sessions are [here](http://huttenhower.sph.harvard.edu/galaxy/u/ignasilucas/h/2020-05-01lefseanalysis)
and [here](http://huttenhower.sph.harvard.edu/galaxy/u/ignasilucas/h/unnamed-history).
Pau suggested to use Lefse to explore the association of amplicon abundances
with binary versions of principal components of life-history data. Some
Acetobacter amplicons (but no genus) was detected significantly associated
mostly with low quality values of PC1 and PC3.

## 2020-04-28
Amparo suggested to include an abundance barplot, that I had avoided before,
because of using too many amplicons, instead of a few taxa. Also, because
the distribution is very skewed, and the plot will not look great. Here I
do some barplots. I also use this folder to make a phylogenetic tree of the
Acetobacter genus, highlighting the position of the amplicons declared
significant by Lefse on `2020-05-01`. Find a preview of the report
[here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-28/figures.html)

## 2020-04-27
Here I complete the analysis of diversity by exploring the relationship
between age or isolines and beta diversity. Curiously, neither phylogenetic
beta diversity of ASVs nor beta diversity of genera decrease with age,
while simple presence-absence beta diversity (not phylogenetically informed)
of ASVs does decrease with age. That is, microbiomes from different flies
resemble more each other late than early in life if we look at ASV composition,
but not if we take into account how similar sequences are.

Then I run again an RDA to correlate ASV abundance and life-history traits.
This time I use only early abundances of ASVs that in `2020-03-17` proved to
require an interaction term between age and isoline to adequately model their
abundances. I apply a logarithmic transformation and keep only ASVs with very
low skewness in their distribution of abundance across isolines (an RDA requirement).
The 18 ASVs left prove significant in a global test, and I proceed with model
selection. There seems to be two Acetobacter clades with negative effects on
either life span or functional aging.

See a preview of the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-27/guts.html).

## 2020-04-25
Summary for group members.

## 2020-04-22
After updating results from the differential abundance analysis (2020-03-17),
which also required updating the phyloseq analysis (2020-02-27), I made sure to
use the same amplicon names in all datasets, and could merge the differential
abundance results with other results, like diversity measures, in order to
have more information about the amplicons. The idea was to improve the RDA
by providing better quality, less noisy information. For example, the diversity
measures, which had not been included in the RDA on 2020-04-21. Here I
compute Faith's phylogenetic diversity index, which turns out to be very
different between early and late samples.

Unfortunately, I have not been able to define, a priori, a set of variables
based on amplicon abundances that can significantly explain variation in
life-history traits through RDA. All models fail the first step, the global
test. I can use the ordistep() function to let the algorithm pick the most
relevant predictor variables, if I relax the significance threshold a bit.
But I don't trust those results. In any case, even doing that, diversity
does not come up as a significant variable.

The report is [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-22/diversity.html)

## 2020-04-21
Full exploration of the RDA results, including all imputed datasets. The
results are not positive. Only one principal component of log-transformed
relative abundances is marginally significantly contributing to variation
in life-history traits. But the amplicons that most heavily load on that
component fail to explain any variation in life-history traits. If the
marginally significant component is not a false positive, we can conclude
that:

* Gut microbiome composition may affect life-history traits, but not
because of the effect of any singular species. The effect, if real,
seems to be combination of taxa.
* The life-history traits most likely to be influenced by gut microbiome
are: the functional measures (climbing speed and its decline with age),
and to some extent reproductive senescence. Neither the life-span nor
the acceleration of mortality rate with age are much affected.

You can see the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-21/interpretation.html)

Some additional work can be done to better characterize the marginally
significant component of amplicon abundances. Maybe some kind of enrichment
analysis could shed some light.

## 2020-04-13
Upon realization that RDA cannot use all the amplicons' abundances, I resort to
reducing the dimensionality of the dataset with PCA, and then use a subset of
the principal components as explanatory variables in RDA. In order to keep the
two isolines with missing values in one of the interesting life-history traits,
I use multiple imputations. Preliminar results suggest there is not much variance
in life-history traits that can be explained by variation in microbiome composition.

See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-13/RDA.html)
About how I chose the logarithmic transformation of relative abundance data, see
[this](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-04-13/transformations.html)


## 2020-03-25
Detailed ordination analysis. Here I learned about ordination analysis. I tested different
distance measures to represent the variation among samples in two dimensions. I confirm
that presence-absence data, rather than quantitative measures of abundance, allow a clearer
distinction of what we consider potentially important sources of variation in two perpendicular
axes: the age of the host flies, and the sequencing batch. I also convinced myself that
double taxon absences should contribute to similarity among communities in our analysis,
despite this not being the case in traditional applications of ordination in ecology.
However, counting or not double absences does not change the display of samples in two
dimensions.

Here I also decide that a redundancy analysis is what we need to determine how and to
what extent bacterial taxa in the gut microbiome affect variation of life-history traits.
I just set up the data matrices and make the analysis run, without getting into the
details. See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-03-25/ordination.html)

## 2020-03-17
I use the DESeq2 package to test for different abundances between early and late
microbiomes. At first, I followed a tutorial for phyloseq users, with very basic
DESeq2 instructions. After reading more about DESeq2, I realize it offers likelihood
ratio tests to compare different models. I have not seen this feature in EdgeR,
for example. It makes me think that the right approach, instead of fitting the
same model to all amplicons, is to separate first the amplicons by their most
adequate model: those that respond to age from those that don't, and so on. Then,
the most interesting amplicons are those that deserve the isoline term in the
abundance model, or even better, those with an interaction term between isoline
and age. Note that applying an LRT is like testing the significance of several
terms simultaneously: I don't care if isoline 39 has an effect, but if any isoline
does. 

See [this](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-03-17/deseq2.html)

## 2020-02-27
I use phyloseq to combine abundance and taxonomy information, and to produce some
plots. See the report [here](https://htmlpreview.github.io/?https://github.com/IgnasiLucas/fly/blob/master/results/2020-02-27/phyloseq.html)


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
