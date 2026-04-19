CC=gcc
CFLAGS=-Wall -O2 -fopenmp
LDFLAGS=-lm -fopenmp

all: galsim

galsim: galsim.c
	$(CC) $(CFLAGS) -o galsim galsim.c $(LDFLAGS)

clean:
	rm -f galsim *.o