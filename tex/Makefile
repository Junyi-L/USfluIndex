cdcflu_Forecast.pdf: references_Forecast.bib

%.pdf: %.tex
	R --vanilla --slave -e "tools::texi2pdf('$<', clean = TRUE)"

%.tex: %.Rnw
	R --slave --no-save --no-restore -e "knitr::knit('$<')"
