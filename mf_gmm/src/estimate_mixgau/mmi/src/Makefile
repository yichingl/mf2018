CC = gcc
#CC = /opt/turing/intel/bin/icc
LIBS = -L/homes/bhiksha/library/$(arch)
INCL = -I/homes/bhiksha/library
LDFLAGS = -lm 
LDFLAGSDBG = -lm
CFLAGS =  -O3 -DIGNORE_READ_ERR -ffast-math

SOURCES = allocs.c  gauss.c  mgau_em.c  mgau_mmi.c vq.c \
          store_distribution.c  twoclassmain.c

TARGET = $(arch)/mgaummi
TARGETDBG = $(arch)/mgaummi.dbg

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LIBS) $(LDFLAGSDBG) -o $(TARGETDBG)

opt: $SOURCES
	$(CC) $(CFLAGS) $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)
	make clean

$SOURCES: header.h

clean:
	rm -f *.o core norm *.c.* *.u
