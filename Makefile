all:
	gcc -o bin/c-anderson `gsl-config --libs --cflags` -lm src/anderson.c

clean:
	rm bin/c-anderson
