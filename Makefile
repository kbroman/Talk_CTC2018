R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

docs/index.html: slides.Rmd docs/kbroman.css R/hs_fig.R R/gve_mixup_scheme.R
	cp slides.Rmd docs/index.Rmd
	cd docs;R $(R_OPTS) -e "rmarkdown::render('index.Rmd')"
	rm docs/index.Rmd

docs/kbroman.css: kbroman.css
	cp kbroman.css docs/

ctc2018.zip: docs/index.html docs/broman_ctc2018.pdf
	mkdir docs/ctc2018
	cd docs;cp -r *.html *.svg *.css index_files/ libs/ ctc2018/
	cp docs/broman_ctc2018.pdf docs/ctc2018/
	cd docs;zip -r ctc2018.zip ctc2018/
	rm -r docs/ctc2018/
	mv docs/ctc2018.zip .

docs/broman_ctc2018.pdf: docs/index.html
	R -e "file <- paste0('file://', normalizePath('docs/index.html'));webshot::webshot(file, '$@')"
	mv broman_ctc2018.pdf docs/

all: ctc2018.zip
zip: ctc2018.zip
pdf: docs/broman_ctc2018.pdf
