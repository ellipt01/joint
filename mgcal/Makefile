# MGCAL
# Magnetic and gravitational field calculation library

include ../make.config

DESTDIR        ?= ..

DESTLIBDIR      = $(DESTDIR)/lib
LIBS            = -lm $(OPENMP_FLG)
CXXFLAGS        = -I./include

LIBSRC_OBJS     = src/magnetic.o src/gravity.o src/grid.o src/kernel.o\
                  src/vector3d.o src/data_array.o src/io.o src/mgcal.o\
                  src/source.o src/private/util.o

DEPS = $(LIBSRC_OBJS:.o=.d)

all:            libmgcal.a

libmgcal.a:     $(LIBSRC_OBJS)
			$(AR) r $@ $(LIBSRC_OBJS)
.c.o:
			$(CC) $(CFLAGS) -MMD -MP -o $*.o -c $(CXXFLAGS) $< $(OPENMP_FLG)

# =========================
# 追加: install
# =========================
install:        libmgcal.a
			@mkdir -p $(DESTLIBDIR)
			$(INSTALL) -m 644 libmgcal.a $(DESTLIBDIR)

clean:
			$(RM) $(DEPS)
			@ for i in $(LIBSRC_OBJS) ; do \
				$(RM) $$i ; \
			done
			$(RM) libmgcal.a
			$(RM) *~ */*~

uninstall:
			$(RM) $(DESTLIBDIR)/libmgcal.a

-include $(DEPS)
