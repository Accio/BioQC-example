## Github does not render svg files in .md documents. We use the rawgit.com service
## to serve the svg files with the correct content-type header and display them
## within the markdown documents
svgfilter=sed -i -r 's/"(.+)\.svg"/"https:\/\/cdn.rawgit.com\/accio\/BioQC-example\/master\/\1\.svg"/g'

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

