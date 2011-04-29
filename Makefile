all: compute/c-anderson
.PHONY: all

compute/c-anderson:
	gcc -o compute/c-anderson `gsl-config --libs --cflags` -lm compute/c-anderson.c

clean:
	rm compute/c-anderson
.PHONY: clean
