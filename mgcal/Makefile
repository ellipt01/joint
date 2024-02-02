# MGCAL
# Magnetic and gravitational field calculation library

include ../make.config

DESTLIBDIR	= ../lib
LIBS		= -lm $(OPENMP_FLG)
CPPFLAGS	= -I./include

LIBSRC_OBJS	= src/magnetic.o src/gravity.o src/grid.o src/kernel.o\
			  src/vector3d.o src/data_array.o src/io.o src/mgcal.o\
			  src/source.o src/private/util.o

all	:		libmgcal

libmgcal:	$(LIBSRC_OBJS)
			$(AR) r $(DESTLIBDIR)/$@.a $(LIBSRC_OBJS)

.c.o:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

install:

clean-objs:
			@ for i in $(LIBSRC_OBJS) ; do \
				$(RM) $$i ; \
			done
			$(RM) *~ */*~

clean:		clean-objs
			$(RM) $(DESTLIBDIR)/libmgcal.a

