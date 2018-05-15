arch=linux
CC = gcc
LDFLAGS = -lgslcblas -lm 
LDFLAGSDBG = -lgslcblas -lm
CFLAGS =  -O3 -ffast-math

SOURCES = extras.c bagofwords.c computebagofwords.c matops.c
TARGET = $(arch)/bow
TARGETDBG = $(arch)/bow.dbg

opt: $SOURCES
	$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(TARGET)

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(SOURCES) $(LDFLAGSDBG) -o $(TARGETDBG)

$SOURCES: header.h

clean:
	rm -f *.o core norm 
