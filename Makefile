include make.config

DESTDIR		= ./bin
DESTLIBDIR	= ./lib

LOCALLIBS	= -L./lib -lmgcal -lmmreal
LIBS		= $(BLAS_LIB) -lm
CFLAGS		= -O3
CPPFLAGS	= -I. -I./include $(BLAS_CFLAGS)\
			  -I./mgcal/include -I./mmreal/include $(OPENMP_FLG)

COMMON_OBJS	= 

JINV_OBJS	= src/main.o src/Joint.o src/Inversion.o src/ADMM.o src/mADMM.o\
			  src/Kernel.o src/utils.o

OBJS		= $(JINV_OBJS) $(COMMON_OBJS)

SUBDIRS		= mgcal mmreal

PROGRAMS	= jinv

all	:		$(SUBDIRS) $(PROGRAMS)

jinv:		$(JINV_OBJS) $(COMMON_OBJS)
			$(CPP) $(CFLAGS) -o $@ $(JINV_OBJS) $(COMMON_OBJS) $(CPPFLAGS) $(LOCALLIBS) $(LIBS)

$(SUBDIRS):	FORCE
			$(MAKE) -C $@

FORCE:

.c:
			$(CC) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

.cc.o:
			$(CPP) $(CFLAGS) -o $*.o -c $(CPPFLAGS) $< $(OPENMP_FLG)

install:
			@ if [ ! -d $(DESTDIR) ]; then \
				mkdir $(DESTDIR) ;\
			fi
			@ for i in $(PROGRAMS) ; do \
				echo $(INSTALL) -m 755 $$i $(DESTDIR) ;\
				$(INSTALL) -m 755 $$i $(DESTDIR) ; \
			done

clean-objs:
			$(RM) $(OBJS)
			@ for i in $(OBJS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(PROGRAMS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				$(MAKE) clean-objs -C $$i ; \
			done
			$(RM) *~ */*~


clean:		clean-objs
			@ for i in $(PROGRAMS) ; do \
				$(RM) $(DESTDIR)/$$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				echo $(MAKE) clean -C $$i ; \
				$(MAKE) clean -C $$i ; \
			done

