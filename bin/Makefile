CXX=g++
# SPECIALFLAGS=--thread_safe --exceptions
SPECIALFLAGS=-O2 -g3
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)
GSLLIBS = -lgsl -lgslcblas -lm
SRT_TOP=.
GSL_DIR=.

# CFLAGS = $(SPECIALFLAGS)  -D__DONT_USE_CAFE__ -I$(GSL_DIR)/include -I$(SRT_TOP)/include
CFLAGS = $(SPECIALFLAGS)  -D__DONT_USE_CAFE__ -I./
LFLAGS = $(SPECIALFLAGS) -L./ -L$(SRT_TOP)/lib/$(SRT_SUBDIR) -L$(GSL_DIR)/lib/ -Ltop_dilepton_madgraph/lib/
# -L../madgraph/ -L/d0usr/products/lhapdf/Linux-2-4/v5_1-GCC3_4_3/lib/

RCXX = $(CFLAGS) $(ROOTCFLAGS)
RLXX = -L/usr/lib/x86_64-linux-gnu/root5.34 $(ROOTLIBS) $(GSLLIBS) -ltop_dilepton_madgraph -lLHAPDF
# -ltop_dilepton_madgraph
# $(CERNLIBS)  -lmatrixelement -ldhelas3 -lLHAPDF 

OBJECTS = base_event.o matrix_calibration.o matrix_element.o matrix_ensemble.o matrix_event.o matrix_kinematic_solver.o matrix_me_dxsec.o matrix_me_sample.o matrix_me_totxsec.o matrix_mwt_sample.o matrix_parameters.o matrix_pdfs.o matrix_resolutions.o matrix_sample.o matrix_weighter.o sigma_function.o ttdilepsolve.o

OBJDIR = obj
OBJECTS_WITH_DIR = $(addprefix $(OBJDIR)/,$(OBJECTS))

VPATH = $(OBJDIR):src:bin

%.o: %.cpp
	$(CXX) $(RCXX) -c $< -o $(OBJDIR)/$(@F)

%.o: %.f
	g77 $(SPECIALFLAGS) -c $< -o $(OBJDIR)/$(@F)

%.o: %.c
	gcc $(SPECIALFLAGS) -c $< -o $(OBJDIR)/$(@F)

all: top_dilepton_me_x

# ${OBJECTS}
top_dilepton_me_x: top_dilepton_me_x.o ctype.o ${OBJECTS}
	$(CXX) $(LFLAGS) -o bin/top_dilepton_me_x obj/top_dilepton_me_x.o  obj/ctype.o ${OBJECTS_WITH_DIR} $(RLXX)

cleanobj:
	rm -f *~ *.o $(OBJDIR)/*.o

clean: cleanobj
	rm -f bin/top_dilepton_me_x
