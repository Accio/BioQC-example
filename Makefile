R=R

all:bioqc-kidney.pdf bioqc-simulation.pdf

bioqc-simulation.pdf:bioqc-simulation.tex
	pdflatex bioqc-simulation.tex
	pdflatex bioqc-simulation.tex

bioqc-simulation.tex:bioqc-simulation.Rnw
	${R} CMD Sweave bioqc-simulation.Rnw

bioqc-kidney.pdf:bioqc-kidney.tex
	pdflatex bioqc-kidney.tex
	pdflatex bioqc-kidney.tex

bioqc-kidney.tex:bioqc-kidney.Rnw
	${R} CMD Sweave bioqc-kidney.Rnw
