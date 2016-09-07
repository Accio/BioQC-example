R=R

all:bioqc-kidney.html bioqc-simulation.html

bioqc-simulation.html:bioqc-simulation.Rmd
	 Rscript -e "rmarkdown::render('bioqc-simulation.Rmd')"

bioqc-kidney.html:bioqc-kidney.Rmd
	Rscript -e "rmarkdown::render('bioqc-kidney.Rmd')"

