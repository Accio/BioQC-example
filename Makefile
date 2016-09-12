R=R

all:bioqc-kidney.html bioqc-simulation.html

clean:
	rm -f bioqc-kidney.html bioqc-simulation.html
	rm -f bioqc-kidney.md bioqc-simulation.md
	rm -rf *_files
	rm -rf *_cache

bioqc-simulation.html:bioqc-simulation.Rmd
	 Rscript -e "rmarkdown::render('bioqc-simulation.Rmd')"

bioqc-kidney.html:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd')"

