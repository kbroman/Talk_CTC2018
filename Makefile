slides.html: slides.Rmd kbroman.css
	cp $^ doc/
	cd doc;R -e "rmarkdown::render('$<')"
	rm doc/$<
