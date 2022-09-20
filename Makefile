timing.md: timing.Rmd
	Rscript -e "rmarkdown::render('timing.Rmd', rmarkdown::md_document(preserve_yaml = TRUE, variant = 'gfm', pandoc_args = '--markdown-headings=atx'))"

efficiency.md: efficiency.Rmd
	Rscript -e "rmarkdown::render('efficiency.Rmd', rmarkdown::md_document(preserve_yaml = TRUE, variant = 'gfm', pandoc_args = '--markdown-headings=atx'))"

## atx headers ensures headers are all like #, ##, etc. Shouldn't be necessary as of pandoc >= 2.11.2
## 'gfm' ensures that the 'r' tag is put on chunks, so code coloring/highlighting will be done when html is produced.

## use quarto as having trouble with the LaTeX equations when rendering via rmarkdown
#efficiency.md: efficiency.Rmd
#	quarto render efficiency.Rmd --to markdown 

tmp.md: tmp.Rmd
	quarto render tmp.Rmd --to markdown 

