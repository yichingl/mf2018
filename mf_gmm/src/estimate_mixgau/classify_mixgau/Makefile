CC = gcc
arch = linux64
#CC = /opt/turing/intel/bin/icc
##LIBS = -L/usr0/bhiksha/library/$(arch)
##INCL = -I/usr0/bhiksha/library
LDFLAGS = -lm 
CFLAGS =  -O3 -ffast-math -lm -Wall

SOURCES = classify.c allocs.c io.c
SOURCESDELTA = classify.delta.c allocs.c io.c

TARGET = $(arch)/classify
TARGETDELTA = $(arch)/classifydelta
TARGETDBG = $(arch)/classify.dbg

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGETDBG)

opt: $(SOURCES)
	$(CC) $(CFLAGS) $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)
	make clean

optdelta: $(SOURCESDELTA)
	$(CC) $(CFLAGS) $(INCL) $(SOURCESDELTA) $(LIBS) $(LDFLAGS) -o $(TARGETDELTA)
	make clean

clean:
	rm -f *.o core norm 
