all: slides.html

%.html: %.Rmd kbroman.css
	R -e "rmarkdown::render('$<')"
