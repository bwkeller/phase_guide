guide.pdf: guide.aux guide.bbl
	latex guide.tex
	dvipdf guide.dvi

guide.aux: guide.tex phase_guide.eps
	latex guide.tex

guide.bbl: guide.aux
	latex guide.tex
	bibtex guide.aux

phase_guide.eps: gen_phase_guide.py
	./gen_phase_guide.py eps

clean:
	rm guide.dvi guide.aux guide.bbl guide.blg guide.log
