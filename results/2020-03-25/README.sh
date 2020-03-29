#!/bin/bash
#
#				2020-03-25
#				----------
#

if [ ! -e ordination.html ]; then
R --no-save -q -e "render_report <- function(){ \
                                rmarkdown::render('ordination.Rmd', \
                                   output_file = 'ordination.html') \
                               }" \
                  -e "render_report()"
fi
