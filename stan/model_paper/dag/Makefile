model.pdf : model.dot
	dot -Txdot model.dot | dot2tex --template template.tex --autosize --crop > model.tex && pdflatex model.tex
