CC = gcc
#CC = /opt/turing/intel/bin/icc
LIBS = -L/homes/bhiksha/library/$(arch)
INCL = -I/homes/bhiksha/library
LDFLAGS = -lm 
LDFLAGSDBG = -lm 
CFLAGS =  -O3 -DIGNORE_READ_ERR -ffast-math
#CFLAGS =  -O3 -DIGNORE_READ_ERR -fast -parallel

SOURCES = EM_diagonal_covariances.c read_cep.c

TARGET = $(arch)/estimate_mgau
TARGETDBG = $(arch)/estimate_mgau.dbg

debug dbg: $(SOURCES)
	gcc -Wall -g -DDEBUG $(INCL) $(SOURCES) $(LIBS) $(LDFLAGSDBG) -o $(TARGETDBG)

opt: $(SOURCES)
	$(CC) $(CFLAGS) $(INCL) $(SOURCES) $(LIBS) $(LDFLAGS) -o $(TARGET)
	make clean

clean:
	rm -f *.o core
