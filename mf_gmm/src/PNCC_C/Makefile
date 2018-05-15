CC      = gcc
CFLAGS  = -Wall -O2 -I.
LDFLAGS = -lm

NBTARGET = wav2pncc.nb
WBTARGET = wav2pncc.wb
SRCS   = main.c pncc.c io.c matrix.c fft.c dct.c hamming.c

wideband: $(SRCS) pncc.h matrix.h sigproc.h Makefile
	$(CC) -DWIDEBAND=1 -o $(WBTARGET) $(CFLAGS) $(SRCS) $(LDFLAGS)

narrowband: $(SRCS) pncc.h matrix.h sigproc.h Makefile
	$(CC) -o $(NBTARGET) $(CFLAGS) $(SRCS) $(LDFLAGS)

clean:
	/bin/rm *.o
