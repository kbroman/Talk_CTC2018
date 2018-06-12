R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

docs/index.html: slides.Rmd kbroman.css
	cp $^ docs/
	mv docs/slides.Rmd docs/index.Rmd
	cd docs;R $(R_OPTS) -e "rmarkdown::render('index.Rmd')"
	rm docs/index.Rmd
