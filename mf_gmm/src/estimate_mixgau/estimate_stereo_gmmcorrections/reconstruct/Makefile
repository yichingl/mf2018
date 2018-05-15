arch = linux64
CC = gcc
LDFLAGS = -lm
CFLAGS =  -O3 

SOURCES = lib.c  IO.c  estimate_corrections.c  main.c

TARGET = $(arch)/em_stereo
TARGETDBG = $(arch)/em_stereo.dbg

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LDFLAGS) -o $(TARGETDBG)

opt: $(SOURCES)
	gcc -Wall -O3 $(INCL) $(SOURCES) $(LDFLAGS) -o $(TARGET)
	make clean

$SOURCES: header.h

clean:
	rm -f *.o core norm *.c.* *.u
