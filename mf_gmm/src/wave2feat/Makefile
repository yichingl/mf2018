VPATH = .

arch = alpha_osf1
#arch = linux
INSTD = /usr0/mseltzer/src/new_fe_wrapper
OBJPATH = /usr0/mseltzer/src/new_fe_wrapper
SRCPATH = /usr0/mseltzer/src/new_fe_wrapper

sun_CFLAGS = -Bstatic -O4 -fsingle -dalign  -DSUN=1
hpux_CFLAGS = -DMEM_DEBUG=0 -DHPUX
alpha_osf1_CFLAGS = -O2 -Dalpha -Dint32=int -Dint16=short \
		-Dfloat32=float -Dfloat64=double
linux_CFLAGS = -m32 -Dintel -Dint32=int -Dint16=short \
		-Dfloat32=float -Dfloat64=double

alpha_osf1_CC = gcc
linux_CC = gcc


CFLAGS = ${$(arch)_CFLAGS} \
	-I$(VPATH) 




# Fri Mar 10 18:06:15 EST 1995 by PJM:
#
# make sure you call the libs in this order (1st sp, 2nd util_nist) , 
# otherwise it fails!!

LIBS = -lm 
TARGET = wave2feat
#test_newfe
OBJS  = 	\
	new_fe_wrapper.o new_fe.o new_fe_sp.o
SRCS  = 	\
	new_fe_wrapper.c new_fe.c new_fe_sp.c
default:
	( $(MAKE) -f $(VPATH)/Makefile CC=${$(arch)_CC}  $(TARGET))


$(TARGET): $(OBJS) 
	$(CC) -o $@ -O2 $(CFLAGS) $(OBJS) $(LIBS)
	strip $(TARGET)
#	mv $(TARGET) $arch/$(TARGET)

debug: $(SRCS) 
	$(CC) -o $(TARGET).$@ -g $(CFLAGS) $(SRCS) $(LIBS)


clean:
	rm -f *.o core  *.c.*

#.c.o:
#	$(CC) -c -o $@ $< $(CFLAGS)

