CFLAGS=-Wall -Ofast -march=native
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: galsim.o graphics.o displaygal.o
	gcc -o galsim galsim.o graphics.o displaygal.o $(LDFLAGS)

galsim.o: galsim.c graphics.h displaygal.h
	gcc $(CFLAGS) $(INCLUDES) -c galsim.c

displaygal.o: displaygal.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c displaygal.c

graphics.o: graphics.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics.c

clean:
	rm -f ./galsim *.o *.gal