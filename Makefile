R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

docs/index.html: slides.Rmd kbroman.css R/hs_fig.R R/gve_mixup_scheme.R
	cp slides.Rmd docs/index.Rmd
	cp kbroman.css docs/
	cd docs;R $(R_OPTS) -e "rmarkdown::render('index.Rmd')"
	rm docs/index.Rmd
