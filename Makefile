all:bioqc-kidney.md bioqc-simulation.md bioqc-wmw-test-performance.md

clean:
	rm -fv bioqc-kidney.md bioqc-simulation.md
	rm -fv bioqc-kidney.md bioqc-simulation.md
	rm -rfv *_files
	rm -rfv *_cache

bioqc-simulation.md:bioqc-simulation.Rmd
	Rscript -e "rmarkdown::render('bioqc-simulation.Rmd', output_format='all')"

bioqc-kidney.md:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd', output_format='all')"

bioqc-wmw-test-performance.md:bioqc-wmw-test-performance.Rmd
	Rscript -e "rmarkdown::render('bioqc-wmw-test-performance.Rmd', output_format='all')"
