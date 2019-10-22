## master tex to pdf
cdcflu_Forecast.pdf: cdcflu_Forecast.tex references_Forecast.bib
	R --vanilla --slave -e "tools::texi2pdf('$<', clean = TRUE)"

## master Rnw to tex
cdcflu_Forecast.tex: cdcflu_Forecast.Rnw ../*.r
	R --slave --no-save --no-restore -e "knitr::knit('$<')"
