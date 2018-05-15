CC      = gcc
CFLAGS  = -Wall -O2 -I.
LDFLAGS = -lm

TARGET = wav2pncc
SRCS   = main.c pncc.c io.c matrix.c fft.c dct.c hamming.c
OBJS   = $(SRCS:.c=.o)

$(TARGET): $(OBJS) Makefile
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LDFLAGS)

.c.o: pncc.h matrix.h sigproc.h Makefile
	$(CC) $(CFLAGS) -c $<

clean:
	-rm $(OBJS)
