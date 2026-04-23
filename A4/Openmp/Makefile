CC=gcc
CFLAGS=-Wall -O3 -fopenmp -fopt-info-vec
LDFLAGS=-lm -fopenmp

all: galsim

galsim: galsim.c
	$(CC) $(CFLAGS) -o galsim galsim.c $(LDFLAGS)

clean:
	rm -f galsim *.o