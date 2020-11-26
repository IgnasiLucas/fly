#!/bin/bash
#
#				2020-04-13
#				==========
#

if [ ! -e transformations.html ]; then
   R -q --no-save -e "rmarkdown::render('transformations.Rmd')"
fi

if [ ! -e RDA.html ]; then
   R -q --no-save -e "rmarkdown::render('RDA.Rmd')"
fi
