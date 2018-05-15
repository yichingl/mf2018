arch=linux
CC = gcc
LDFLAGS = -lgslcblas -lm 
LDFLAGSDBG = -lgslcblas -lm
CFLAGS =  -O3 -ffast-math

SOURCES = allocs.c gaussrand.c io.c main.c matops.c vq.c 
TARGET = $(arch)/vq
TARGETDBG = $(arch)/vq.dbg

opt: $SOURCES
	$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(TARGET)

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(SOURCES) $(LDFLAGSDBG) -o $(TARGETDBG)

$SOURCES: header.h

clean:
	rm -f *.o core norm 
