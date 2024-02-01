include ../make.config

DESTLIBDIR	= ../lib
LIBS		= -lm $(OPENMP_FLG)
CPPFLAGS	= -I./include

LIBSRC_OBJS	= src/mmio.o src/mmreal.o

all	:		libmmreal

libmmreal:	$(LIBSRC_OBJS)
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
			$(RM) $(DESTLIBDIR)/libmmreal.a

