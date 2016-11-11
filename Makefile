
R=R

all: setup bioqc-simulation.md bioqc-kidney.md bioqc-wmw-test-performance.md create_gtex_signatures.md compare_signatures.md 

setup:
	mkdir -p gmt

clean:
	rm -fv *.html
	find . -maxdepth 1 -type f -name "*.md" | grep -v README | xargs rm -fv
	rm -rfv *_files
	rm -rfv fig/*
	rm -rfv gmt/*

wipe: clean
	rm -rfv *_cache

create_gtex_signatures.md: create_gtex_signatures.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

compare_signatures.md: compare_signatures.Rmd
	mkdir -p compare_signatures_files
	Rscript -e "rmarkdown::render('$<', output_format='all')"

bioqc-simulation.md:bioqc-simulation.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

bioqc-kidney.md:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"

bioqc-wmw-test-performance.md:bioqc-wmw-test-performance.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='all')"
