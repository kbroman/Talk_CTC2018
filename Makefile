slides.html: slides.Rmd kbroman.css
	cp $^ docs/
	cd docs;R -e "rmarkdown::render('$<')"
	rm docs/$<
