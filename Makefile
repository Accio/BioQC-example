
R=R

all: setup create_gtex_signatures.html compare_signatures.html

setup:
	mkdir -p gmt cache_data fig

clean:
	rm -fv *.html
	rm -rfv *_files
	rm -rfv fig/*
	rm -rfv gmt/*

wipe: clean
	rm -rfv *_cache
	rm -rfv cache_data/*

create_gtex_signatures.html: create_gtex_signatures.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

compare_signatures.html: compare_signatures.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

bioqc-simulation.md:bioqc-simulation.Rmd
	Rscript -e "rmarkdown::render('bioqc-simulation.Rmd', output_format='all')"

bioqc-kidney.md:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd', output_format='all')"

bioqc-wmw-test-performance.md:bioqc-wmw-test-performance.Rmd
	Rscript -e "rmarkdown::render('bioqc-wmw-test-performance.Rmd', output_format='all')"
