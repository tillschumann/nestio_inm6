# Makefile for latex documents

MFILE=report

all:
	pdflatex $(MFILE)
	bibtex $(MFILE)
	pdflatex $(MFILE)
	pdflatex $(MFILE)

clean:
	rm -rf *~ *.log *.aux *.bbl *.blg *.idx *.ilg *.ind *.lof *.out *.pfg *.toc

mrproper: clean
	rm -rf *.pdf *.ps *.dvi 
