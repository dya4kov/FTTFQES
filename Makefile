IDIR=./hdr
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR)
ODIR=./obj
BINDIR=./bin
SDIR=./src
LIBS=-lgsl -lgslcblas -lm

_DEPS=spec_func.h
      
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=main.o \
     spec_func.o

OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

all: makeobjdir makebindir $(OBJ) TFCS

makeobjdir:
	mkdir $(ODIR)

makebindir:
	mkdir $(BINDIR)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

TFCS: $(OBJ) 
	$(CC) -o $(BINDIR)/$@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -rf $(ODIR) *~ $(IDIR)/*~ $(SDIR)/*~ $(BINDIR)
