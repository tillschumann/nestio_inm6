# Makefile for latex documents

MFILE=schumann_yic_2015_4paper.tex

all:
	pdflatex $(MFILE)
	pdflatex $(MFILE)
	pdflatex $(MFILE)

clean:
	rm -rf *~ *.log *.aux *.bbl *.blg *.idx *.ilg *.ind *.lof *.out *.pfg *.toc

mrproper: clean
	rm -rf *.pdf *.ps *.dvi 
