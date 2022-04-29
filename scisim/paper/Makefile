name=paper

all: 
	pdflatex $(name).tex
	bibtex $(name)
	pdflatex $(name).tex
	pdflatex $(name).tex

.PHONY : clean
clean :
	rm $(name).aux
	rm $(name).log
	rm $(name).out
	rm $(name).blg
	rm -f $(name).fls
	rm -f $(name).fls_latexmk

