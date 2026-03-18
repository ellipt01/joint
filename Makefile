include make.config

DESTDIR         ?= .

LOCALLIBS       = -L./mgcal -lmgcal
LIBS            = $(BLAS_LIB) -lm
CFLAGS          = -O3
CXXFLAGS        = -I. -I./include $(BLAS_CFLAGS)\
			  -I./mgcal/include

COMMON_OBJS	= 

JINV_OBJS       = src/main.o src/Joint.o src/ADMM.o src/mADMM.o src/Kernel.o

OBJS            = $(JINV_OBJS) $(COMMON_OBJS)

SUBDIRS		= mgcal

PROGRAMS        = jinv

DEPS = $(OBJS:.o=.d)

all	:		$(SUBDIRS) $(PROGRAMS)

jinv:			$(JINV_OBJS) $(COMMON_OBJS)
			$(CXX) $(CFLAGS) -o $@ $(JINV_OBJS) $(COMMON_OBJS) $(CXXFLAGS) $(LOCALLIBS) $(LIBS) $(OPENMP_FLG)

$(SUBDIRS):	FORCE
			$(MAKE) -C $@

FORCE:

.c.o:
			$(CC) $(CFLAGS) -MMD -MP -o $*.o -c $(CXXFLAGS) $< $(OPENMP_FLG)

.cc.o:
			$(CXX) $(CFLAGS) -MMD -MP -o $*.o -c $(CXXFLAGS) $< $(OPENMP_FLG)

install:
			@mkdir -p $(DESTDIR)/bin
			@ for i in $(PROGRAMS) ; do \
				echo $(INSTALL) -m 755 $$i $(DESTDIR)/bin ;\
				$(INSTALL) -m 755 $$i $(DESTDIR)/bin ; \
			done
			@ for i in $(SUBDIRS) ; do \
				echo $(MAKE) install -C $$i ; \
				$(MAKE) install -C $$i ; \
			done
clean:
			$(RM) $(OBJS)
			$(RM) $(DEPS)
			@ for i in $(OBJS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(PROGRAMS) ; do \
				$(RM) $$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				$(MAKE) clean -C $$i ; \
			done
			$(RM) *~ */*~

uninstall:
			@ for i in $(PROGRAMS) ; do \
				$(RM) $(DESTDIR)/bin/$$i ; \
			done
			@ for i in $(SUBDIRS) ; do \
				echo $(MAKE) uninstall -C $$i ; \
				$(MAKE) uninstall -C $$i ; \
			done

-include $(DEPS)
