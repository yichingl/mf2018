arch=linux
CC = gcc
LDFLAGS = -lgslcblas -lm 
LDFLAGSDBG = -lgslcblas -lm 
CFLAGS =  -O3

SOURCES = allocs.c io.c main.c
TARGET = $(arch)/getcov
TARGETDBG = $(arch)/getcov.dbg

opt: $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(TARGET)

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(SOURCES) $(LDFLAGSDBG) -o $(TARGETDBG)

$SOURCES: header.h

clean:
	rm -f *.o core norm 
