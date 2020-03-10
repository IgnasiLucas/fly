#!/bin/bash
#
#				2020-02-27
#				----------
#
# In 2020-02-11 I created dataframe TA, to combine information on taxonomy and
# abundance in one table. But then I realized that I was re-inventing the wheel.
# Of course there are packages to process both pieces of information. Here I use
# `phyloseq`, following some well documented examples and tutorials.

# SET UP
# ======

DATA='../2020-02-11/taxonomy.RData'
RUN1='../../data/RUNGEN52_2019/fastq'
RUN2='../../data/RUNGEN66/fastq'

# ******************************************************************************

if [ ! -e phyloseq.html ]; then
R --no-save -q -e "render_report <- function(data, run1, run2){ \
                                rmarkdown::render('phyloseq.Rmd', \
                                   params = list(DATA = data, RUN1 = run1, RUN2 = run2), \
                                   output_file = 'phyloseq.html') \
                               }" \
                  -e "render_report('$DATA', '$RUN1', '$RUN2')"
fi
