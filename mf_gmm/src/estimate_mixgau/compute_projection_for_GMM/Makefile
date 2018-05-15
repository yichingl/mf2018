CC = gcc
#CC = /opt/turing/intel/bin/icc
LIBS = -L/homes/bhiksha/library/$(arch)
INCL = -I/homes/bhiksha/library
LDFLAGS = -lm -layesha
CFLAGS =  -O3 -DIGNORE_READ_ERR -ffast-math
#CFLAGS =  -O3 -DIGNORE_READ_ERR -fast -parallel

SOURCES = main.c vq.c em_full_vars.c gauss.c \
          read_backup_distributions.c store_distribution.c 

TARGET = $(arch)/estimate_mgau

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)

opt: $SOURCES
	$(CC) $(CFLAGS) $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)
	make clean

$SOURCES: header.h

clean:
	rm -f *.o core norm *.c.* *.u
