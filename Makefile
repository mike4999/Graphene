all: compute/c-anderson report/report.pdf
.PHONY: all

compute/c-anderson: compute/c-anderson.c
	gcc -o compute/c-anderson `gsl-config --libs --cflags` -lm compute/c-anderson.c

report/report.pdf: report/report.tex
	cd report; pdflatex report.tex; pdflatex report.tex

clean:
	rm -f compute/c-anderson
	rm -f report/report.aux report/report.log report/report.pdf
.PHONY: clean
