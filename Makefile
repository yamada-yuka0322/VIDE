# Compiler choice:

# Gcc
CC = icc
CFLAGS = -O3 -ggdb

MLIBS	=   -lm

#QLIB = -L../qhull-2003.1/src/
#QINC = -I../qhull-2003.1/src/
QLIB = -L../qhull2002.1/src/
QINC = -I../qhull2002.1/src/

###############

all: vozinit voz1b1 voztie jozov

jozov: jozov.o findrtop.o Makefile
	$(CC) $(CFLAGS) -o jozov jozov.o findrtop.o $(MLIBS)
jozov.o: jozov.c Makefile
	$(CC) $(CFLAGS) -c jozov.c
findrtop.o: findrtop.c Makefile
	$(CC) $(CFLAGS) -c findrtop.c

voz1b1: voz1b1.o readfiles.o vozutil.o voz.h Makefile
	$(CC)  -o voz1b1 $(CFLAGS) voz1b1.o readfiles.o vozutil.o -L. $(QLIB) -lqhull $(MLIBS) 
voz1b1.o: voz1b1.c Makefile
	$(CC) $(CFLAGS) $(QINC) -c voz1b1.c
vozutil.o: vozutil.c Makefile
	$(CC) $(CFLAGS) $(QINC) -c vozutil.c

vozinit: vozinit.o readfiles.o voz.h Makefile
	$(CC)  -o vozinit $(CFLAGS) vozinit.o readfiles.o -L. $(MLIBS) 
vozinit.o: vozinit.c Makefile
	$(CC) $(CFLAGS) -c vozinit.c $(QINC)

voztie: voztie.o readfiles.o Makefile
	$(CC)  -o voztie $(CFLAGS) voztie.o readfiles.o
voztie.o: voztie.c Makefile
	$(CC) $(CFLAGS) -c voztie.c

