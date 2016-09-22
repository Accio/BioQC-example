all:bioqc-kidney.html bioqc-simulation.html

clean:
	rm -fv bioqc-kidney.html bioqc-simulation.html
	rm -fv bioqc-kidney.md bioqc-simulation.md
	rm -rfv *_files
	rm -rfv *_cache

bioqc-simulation.html:bioqc-simulation.Rmd
	Rscript -e "rmarkdown::render('bioqc-simulation.Rmd', output_format='all')"

bioqc-kidney.html:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd', output_format='all')"

