#!/bin/bash
#
#       2020-02-11
#       ----------
#
# Here I will assign taxonomy to every amplicon sequence variant.
# I will follow the dada2 tutorial. See details in accompanying
# taxonomy.Rmd file.

TABSEQDATA='../2019-12-12/denoising.RData'
SILVA_TRAINING='../../data/silva_nr99_v138_train_set.fa.gz'
SILVA_SPECIES='../../data/silva_species_assignment_v138.fa.gz'
NUM_THREADS=12

if [ ! -e $SILVA_TRAINING ]; then
   echo "Error: $SILVA_TRINING not found"
   exit
fi

if [ ! -e $SILVA_SPECIES ]; then
   echo "Error: $SILVA_SPECIES not found"
   exit
fi

if [ ! -e taxonomy.html ]; then
   # Below, instead of executing rmarkdown::render directly, I create a function,
   # that takes as arguments the parameters that I need to pass to rmarkdown::render.
   # This is the best way I've found to render an Rmd with parameters. "Why do I pass
   # some values as parameters?" The paths to the database files, for example, are
   # system-dependent. I prefer to set them in only one place, where they are easy
   # to change. Note that the parameter values passed this way will override the
   # values specified in the YAML section of the Rmd file.
   R --no-save -q -e "render_report <- function(tabseqdata, silva_training, silva_species, num_threads){ \
                                rmarkdown::render('taxonomy.Rmd', \
                                   params = list(TABSEQDATA = tabseqdata, \
                                                 SILVA_TRAINING = silva_training, \
                                                 SILVA_SPECIES  = silva_species, \
                                                 NUM_THREADS = num_threads), \
                                   output_file = 'taxonomy.html') \
                               }" \
                  -e "render_report('$TABSEQDATA', '$SILVA_TRAINING', '$SILVA_SPECIES', '$NUM_THREADS')"
fi
