svgfilter=sed -i -r 's/"(.+)\.svg"/"https:\/\/cdn.rawgit.com\/grst\/BioQC-example\/master\/\1\.svg"/g'

all:bioqc-kidney.html bioqc-simulation.html

clean:
	rm -fv bioqc-kidney.html bioqc-simulation.html
	rm -fv bioqc-kidney.md bioqc-simulation.md
	rm -rfv *_files
	rm -rfv *_cache

bioqc-simulation.html:bioqc-simulation.Rmd
	 Rscript -e "rmarkdown::render('bioqc-simulation.Rmd')"
	 $(svgfilter) bioqc-simulation.md

bioqc-kidney.html:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd')"
	$(svgfilter) bioqc-kidney.md

