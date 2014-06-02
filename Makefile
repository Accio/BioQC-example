all:bioqc-kidney.pdf

bioqc-kidney.pdf:bioqc-kidney.tex
	pdflatex bioqc-kidney.tex
	pdflatex bioqc-kidney.tex

bioqc-kidney.tex:bioqc-kidney.Rnw
	R-devel CMD Sweave bioqc-kidney.Rnw