#!/bin/bash
#
# Note on reproducibility
# =======================
#
# I tried to use an Rmd file as the source file, because R can also execute bash code.
# Rendering an Rmd file would produce a more readable version of the source code, at the
# same time that the results are re-generated. But there is a problem: when using rstudio-
# server, all commands run in a default environment that I cannot change from within that
# session. That is, I cannot activate the conda environment required by an R session from
# within that same session. I recall it is also impossible for an R script to render
# itself.
#
# Actually, I cannot run "conda activate <my_env>" from a script. See https://github.com/conda/conda/issues/7980
# The solution here is to use "source activate", instead.
#
# Following a [recommendation](https://kaust-vislab.github.io/introduction-to-conda-for-data-scientists/03-sharing-environments/index.html) 
# by David R. Pugh, I should create a yml file with only the top packages of the dependency
# graph, and install the environment in the current directory. However, qiime2 is installed
# in conda from a different yml depending on the operating system. And installing only
# the package called qiime2 (from channel qiime2) does not work. Plus, I already started
# using qiime on 2019-11-14, and re-creating the environment in every folder would be redundant.
# I use a system-wide enviornment called qiime2-2019.10, created following the instructions
# in the [qiime website](https://docs.qiime2.org/2019.10/install/native/#install-qiime-2-within-a-conda-environment).
# I have to create the yml file with `conda env export`, and it is only expected to work in
# Linux. Sorry. It is saved here as `environment.yml`. You must have activated the qiime2-2019.10
# environment before rendering the source `.Rmd` file.

DATADIR=../2019-11-14/FromMerged
SAMPLE=(E10A E10B E10C       E11A E11B E11C       E12A E12B E12C       E14A E14B      E14D
        E15A      E15C E15D  E17A E17B E17C E17D  E19A E19B E19C E19D  E20A E20B E20C
        E22A E22B      E22D       E23B E23C E23D  E24A E24B E24C       E6A  E6B       E6D
        L10A L10B L10C L10D  L11A L11B L11C L11D  L12A L12B L12C L12D  L14A L14B L14C L14D
        L15A L15B L15C       L17A L17B L17C L17D  L19A L19B L19C L19D  L20A L20B L20CD
        L22A L22B L22C L22D  L23A           L23D  L24A L24B            L6A  L6B  L6C  L6D)

if [ $CONDA_DEFAULT_ENV != qiime2-2019.10 ]; then
   source activate qiime2-2019.10
fi

# We need a metadata file with sample attributes. I assume, from the names of the samples,
# that the number within the name refers to the fly isoline; the first letter indicate if
# sampling happened "Early" or "Late" in life.
if [ ! -e metadata.tsv ]; then
   echo -e "sample-id\tAge\tIsoline"              > metadata.tsv
   echo -e "#q2:types\tcategorical\tcategorical" >> metadata.tsv
   for i in ${SAMPLE[@]}; do
      AGE=${i:0:1}                               # first letter of sample name.
      ISO=$(echo $i | grep -Eo "[[:digit:]]+")   # Also: $(grep -Eo "[[:digit:]]+" <<< $i)
      echo -e "$i\t$AGE\t$ISO"                   >> metadata.tsv
   done
fi

# This creates a visualization report, to be viewd either with qiime, on a local computer,
# or probably from https://view.qiime2.org/ if your browser supports it.
if [ ! -e FeatureTable.qzv ]; then
   qiime feature-table summarize --i-table $DATADIR/FeatureTable.qza \
                                 --o-visualization  FeatureTable.qzv \
                                 --m-sample-metadata-file metadata.tsv
fi

if [ ! -e FeatureTable.biom ]; then
   qiime tools export --input-path $DATADIR/FeatureTable.qza \
                      --output-path .
   mv feature-table.biom FeatureTable.biom
fi

if [ ! -e ExploreData.html ]; then
   R --save -q -e "rmarkdown::render('ExploreData.Rmd', output_file='ExploreData.html')"
fi

