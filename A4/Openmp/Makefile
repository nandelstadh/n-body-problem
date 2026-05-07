CC=gcc
CFLAGS=-Wall -O3 -fopenmp -fopt-info-vec
INCLUDES=-I/opt/X11/include
LDFLAGS=-lm -fopenmp -L/opt/X11/lib -lX11

all: galsim

galsim: galsim.o graphics.o displaygal.o
	gcc -o galsim galsim.o graphics.o displaygal.o $(LDFLAGS)

galsim.o: ../../galsim.c ../../graphics.h ../../displaygal.h
	gcc $(CFLAGS) $(INCLUDES) -c galsim.c

displaygal.o: ../../displaygal.c ../../graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c ../../displaygal.c

graphics.o: ../../graphics.c ../../graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c ../../graphics.c

clean:
	rm -f ./galsim *.o *.gal