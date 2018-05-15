CC = gcc
LIBS = -L/homes/bhiksha/library/$(arch)
INCL = -I/homes/bhiksha/library
LDFLAGS = -lm -layesha
CFLAGS =  -O3 -DIGNORE_READ_ERR 

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
