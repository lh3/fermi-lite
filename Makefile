CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o misc.o \
			bseq.o htab.o bfc.o \
			rle.o rope.o mrope.o rld0.o \
			unitig.o mag.o bubble.o ksw.o
PROG=		fml-test
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fml-test:$(OBJS) test.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bfc.o: htab.h kmer.h internal.h fml.h kvec.h ksort.h
bseq.o: fml.h kseq.h
htab.o: htab.h kmer.h khash.h
misc.o: internal.h fml.h kstring.h rle.h mrope.h rope.h rld0.h
mrope.o: mrope.h rope.h
rld0.o: rld0.h
rle.o: rle.h
rope.o: rle.h rope.h
test.o: fml.h htab.h kmer.h
