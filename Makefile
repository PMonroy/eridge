CC = gcc
CFLAGS= -O2 -Wall -mcmodel=medium
LIBS = -lteem -lm
RM = rm -rf

PARAMS=parameters.in

all: eridge mkNrrd clean

eridge: eridge.o
	$(CC) $(CFLAGS) -o eridge eridge.o $(LIBS)

mkNrrd: mkNrrd.o
	$(CC) $(CFLAGS) -o mkNrrd mkNrrd.o $(LIBS)

eridge.o: eridge.c
	$(CC) $(CFLAGS) -c eridge.c

mkNrrd.o: mkNrrd.c
	$(CC) $(CFLAGS) -c mkNrrd.c

clean:
	$(RM) *.o 
