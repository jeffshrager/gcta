name=model_paper

all: dag/model.pdf gcta_clean.bib $(name).pdf 

$(name).pdf : $(name).tex
	pdflatex $(name).tex
	bibtex $(name)
	pdflatex $(name).tex
	pdflatex $(name).tex

dag/model.pdf : dag/model.dot
	$(MAKE) -C dag

gcta_clean.bib : gcta.bib
	python preprocbib.py

.PHONY : clean

clean :
	rm $(name).aux
	rm $(name).log
	rm $(name).out
	rm $(name).bbl
	rm $(name).blg
