R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

slides.html: slides.Rmd kbroman.css
	cp $^ docs/
	cd docs;R $(R_OPTS) -e "rmarkdown::render('$<')"
	rm docs/$<
