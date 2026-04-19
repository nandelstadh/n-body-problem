CC=gcc
CFLAGS=-Wall -O2
LDFLAGS=-lpthread -lm

all: galsim

galsim: galsim.c
	$(CC) $(CFLAGS) -o galsim galsim.c $(LDFLAGS)

clean:
	rm -f galsim *.o