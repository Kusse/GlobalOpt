CC	= 	cc 
LIBS	=	-lm
CFLAGS	=	-fast
CFLAGS	=	-g
CFLAGS	=	-O3


OBJS	= 	CNmain.o CNlist.o CNio.o CNcna.o CNpairStat.o
cnA.exc:	$(OBJS)
		$(CC) $(CFLAGS) -o cnA.exc $(OBJS) $(LIBS)

CNmain.o: CNmain.c CNmain.h CNpairStat.h
	${CC} -c $(CFLAGS) CNmain.c

CNlist.o: CNlist.c CNmain.h CNmacros.h CNpbc.h
	${CC} -c $(CFLAGS) CNlist.c

CNio.o: CNio.c CNmain.h CNmacros.h
	${CC} -c $(CFLAGS) CNio.c

CNcna.o: CNcna.c CNmain.h
	${CC} -c $(CFLAGS) CNcna.c

CNpairStat.o: CNpairStat.c CNmain.h CNpairStat.h CNmacros.h
	${CC} -c $(CFLAGS) CNpairStat.c
