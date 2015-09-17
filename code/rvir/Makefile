#MRS Makefile template
include Makefile.local

.PHONY : all clean oclean

#C implicit rule: $(CC) -c $(CPPFLAGS) $(CFLAGS)
CPPFLAGS = 
CFLAGS = -Wall

#linking implicit rule: $(CC) $(LDFLAGS) [objects] $(LOADLIBES) $(LDLIBS)
LDFLAGS = 
LOADLIBES = 
LDLIBS = -lm

ALLEXECS = rvir

all: $(ALLEXECS)

rvir: rvir.o parse.o

ALLOBJS = rvir.o parse.o
ALLINCL = rvir.h
ALLSRC := $(ALLOBJS:.o=.c)

rvir.tar.gz: Makefile $(ALLSRC) $(ALLINCL)
	tar zcvf $@ $^

clean:
	$(RM) $(ALLOBJS) $(ALLEXECS)

oclean:
	$(RM) $(ALLOBJS)
