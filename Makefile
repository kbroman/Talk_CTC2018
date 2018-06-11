all: slides.html

%.html: %.Rmd
	R -e "rmarkdown::render('$<')"
