arch=linux
CC = gcc
LDFLAGS = -lgslcblas -lm 
LDFLAGSDBG = -lgslcblas -lm 
CFLAGS =  -O3

SOURCES = lib.c main.c vq.c fileio.c matops.c
TARGET = $(arch)/vq
TARGETDBG = $(arch)/vq.dbg

opt: $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(TARGET)

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(SOURCES) $(LDFLAGSDBG) -o $(TARGETDBG)

simd: $(SOURCES)
	$(CC) -DSIMD=1 $(CFLAGS) $(SOURCES) $(LDFLAGS) -o $(TARGET).simd

simddbg: $(SOURCES)
	$(CC) -g -DSIMD=1 $(SOURCES) $(LDFLAGS) -o $(TARGETDBG).simd

multithread: $(SOURCES)
	$(CC) -DMULTITHREAD=1 $(CFLAGS) $(SOURCES) $(LDFLAGS) -lpthread -o $(TARGET).multithread

multithreaddbg: $(SOURCES)
	$(CC) -DMULTITHREAD=1 -g $(SOURCES) $(LDFLAGS) -lpthread -o $(TARGETDBG).multithread

multithreadsimd: $(SOURCES)
	$(CC) -DMULTITHREAD=1 -DSIMD=1 $(CFLAGS) $(SOURCES) $(LDFLAGS) -lpthread -o $(TARGET).multithread.simd

multithreadsimddbg: $(SOURCES)
	$(CC) -DMULTITHREAD=1 -DSIMD=1 -g $(SOURCES) $(LDFLAGS) -lpthread -o $(TARGETDBG).multithread.simd

$SOURCES: header.h

clean:
	rm -f *.o core norm 
