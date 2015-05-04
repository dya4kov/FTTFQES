IDIR=./hdr
LDIR=./lib
ODIR=./obj
BDIR=./bin
SDIR=./src
CC=g++
CFLAGS=-I$(IDIR)
LIBS=-lm

TF_POT_DEPS=$(IDIR)/FTTFpotential.h \
	  		$(IDIR)/ThermodynamicFunction.h 
TFQE_POT_DEPS=
TF_MOD_DEPS=
TFQE_MOD_DEPS=
TFS_MOD_DEPS=
TF_EL_STATE_DEPS=
TFQES_DEPS=
Y_DEPS=$(IDIR)/Yfunction.h \

TF_POT_OBJ=$(ODIR)/FTTFpotential.o
TFQE_POT_OBJ=$(ODIR)/FTTFQEpotential.o
TF_MOD_OBJ=$(ODIR)/FTTFmodel.o
TFQE_MOD_OBJ=$(ODIR)/FTTFQEmodel.o
TFS_MOD_OBJ=$(ODIR)/FTTFSmodel.o
TF_EL_STATE_OBJ=$(ODIR)/FTTFelectronicStates.o
TFQES_OBJ=$(ODIR)/FTTFQESmodel.o
Y_OBJ=$(ODIR)/Yfunction.o
PRINT_OBJ=$(ODIR)/printer.o
TIMER_OBJ=$(ODIR)/Timer.o

objdir:
	mkdir -p $(ODIR) 

bindir:
	mkdir -p $(BDIR) 

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< -I$(IDIR)

FTTFpot: objdir bindir $(TF_POT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TF_POT_OBJ) $(LIBS)

FTTFpot: DEPS=$(TF_POT_DEPS)

FTTFQEpot: objdir bindir $(TFQE_POT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TFQE_POT_OBJ) $(LIBS)

FTTFQEpot: DEPS=$(TFQE_POT_DEPS)	

FTTFmodel: objdir bindir $(TF_MOD_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TF_MOD_OBJ) $(LIBS)

FTTFmodel: DEPS=$(TF_MOD_DEPS)

FTTFQEmodel: objdir bindir $(TFQE_MOD_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TFQE_MOD_OBJ) $(LIBS)

FTTFQEmodel: DEPS=$(TFQE_MOD_DEPS)

FTTFeStates: objdir bindir $(TF_EL_STATE_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TF_EL_STATE_OBJ) $(LIBS)

FTTFeStates: DEPS=$(TF_EL_STATE_OBJ)

FTTFQES: objdir bindir $(TFQES_MOD_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TFQES_MOD_OBJ) $(LIBS)

FTTFQES: DEPS=$(TFQES_DEPS)

Ytest: objdir bindir $(Y_OBJ) 
	$(CC) -o $(BDIR)/$@ $(Y_OBJ) $(LIBS)

Ytest: DEPS=$(Y_DEPS)

printTest: objdir bindir $(PRINT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(PRINT_OBJ) $(LIBS)

timer: objdir bindir $(TIMER_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TIMER_OBJ) $(LIBS)

.PHONY: clean

clean:
	rm -rf $(ODIR) *~ $(IDIR)/*~ $(SDIR)/*~ $(BDIR)
