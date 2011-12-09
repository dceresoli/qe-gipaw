TEX = user-manual.tex
PDF = $(TEX:.tex=.pdf)

all: $(PDF)

pdf: $(PDF)

$(PDF): $(TEX)
	pdflatex $(TEX)
	pdflatex $(TEX)

clean:
	-rm -f *~ *.log *.dvi *.aux *.end *.out 2>/dev/null

veryclean: clean
	-rm -f $(PDF)

