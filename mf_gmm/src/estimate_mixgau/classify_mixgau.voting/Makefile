CC = gcc
#CC = /opt/turing/intel/bin/icc
LIBS = -L/home/bhiksha/library/$(arch)
INCL = -I/home/bhiksha/library
LDFLAGS = -lm -layesha
LDFLAGSDBG = -lm -layesha_dbg
CFLAGS =  -O3 -ffast-math -layesha -lm

SOURCES = classify.c

TARGET = $(arch)/classify
TARGETDBG = $(arch)/classify.dbg

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LIBS) $(LDFLAGSDBG) -o $(TARGETDBG)

opt: $(SOURCES)
	$(CC) $(CFLAGS) $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)
	make clean

clean:
	rm -f *.o core norm 
