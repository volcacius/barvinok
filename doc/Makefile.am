export TEXINPUTS := $(srcdir):$(TEXINPUTS)
export BIBINPUTS := $(srcdir):$(BIBINPUTS)
export BSTINPUTS := $(srcdir):$(BSTINPUTS)

barvinok.dvi: FORCE
	(cd ..; echo "@GIT_HEAD_VERSION@") > version.tex
	latex barvinok.tex
	bibtex barvinok
	makeindex barvinok.idx
	makeindex barvinok.glo -o barvinok.gls -s nomencl.ist
	glosstex barvinok.aux $(srcdir)/barvinok.gdf
	makeindex barvinok.gxs -o barvinok.glx -s glosstex.ist
	latex barvinok.tex
	latex barvinok.tex

barvinok.pdf: barvinok.dvi
	dvips -Ppdf barvinok.dvi
	ps2pdf barvinok.ps

FORCE:
