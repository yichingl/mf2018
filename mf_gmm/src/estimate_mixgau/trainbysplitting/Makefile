CC = gcc
LDFLAGS = -lm 
CFLAGS =  -D "REAL=float" -D "MAXITER=15"

TARGET = $(arch)/estimate_mgau

debug dbg: 
	gcc -Wall -g -DDEBUG $(CFLAGS) $(INCL) main.c EM_diagonal_covariances.c allocs.c io.c $(LDFLAGS) -o linux/estimate_mgau.dbg -DSKIPBADVECTORS

opt:
	gcc -O2 $(CFLAGS) $(INCL) main.c EM_diagonal_covariances.c allocs.c io.c $(LIBS) $(LDFLAGS) -o linux/estimate_mgau -DSKIPBADVECTORS

debugshared: 
	gcc -Wall -g -DDEBUG $(CFLAGS) $(INCL) main.c EM_sharedcovariance.c allocs.c io.c $(LIBS) $(LDFLAGS) -o linux/estimate_mgaushared.dbg

optshared: 
	gcc -O2 $(CFLAGS) $(INCL) main.c EM_sharedcovariance.c allocs.c io.c $(LIBS) $(LDFLAGS) -o linux/estimate_mgaushared

clean:
	rm -f *.o core norm *.c.* *.u
