#MRS Makefile template
include ./Makefile.local

.PHONY : all clean oclean

#Fortran implicit rule: $(FC) -c $(FFLAGS)
FFLAGS = -Wall 
FLIBS = 

#Fortran linking example:
#target: $(OBJS)
#	$(FC) $(FFLAGS) $^ -o $@ $(FLIBS)

#C implicit rule: $(CC) -c $(CPPFLAGS) $(CFLAGS)
CPPFLAGS = 
CFLAGS = -g -Wall  $(FAST) $(GSL_FLAGS) -lm

#linking implicit rule: $(CC) $(LDFLAGS) [objects] $(LOADLIBES) $(LDLIBS)
LDFLAGS = $(FAST)
LOADLIBES = 
LDLIBS = 

ALLEXECS = stmf

all: $(ALLEXECS)

stmf: sigma2.o mrseishu.o mrsromb.o massfunc.o growth.o bias.o

EISHUOBJS = mrseishu.o

$(EISHUOBJS): eishu.h

PKOBJS = sigma2.o massfunc.o stmf.o

$(PKOBJS): pk.h bessel.h

GROWOBJS = growth.o stmf.o massfunc.o bias.o

$(GROWOBJS): growth.h

ALLOBJS = $(EISHUOBJS) $(PKOBJS) $(GROWOBJS) mrsromb.o
ALLINCL = eishu.h pk.h bessel.h growth.h
ALLSRC := $(ALLOBJS:.o=.c)

stmf.tar: Makefile $(ALLSRC) $(ALLINCL) mac.sm
	tar cvf $@ $^

clean:
	$(RM) $(ALLOBJS) $(ALLEXECS)

oclean:
	$(RM) $(ALLOBJS)
