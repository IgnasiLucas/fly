#!/bin/bash
#
#				2020-04-21
#				==========
#

if [ ! -e interpretation.html ]; then
   R -q --no-save -e "rmarkdown::render('interpretation.Rmd')"
fi
