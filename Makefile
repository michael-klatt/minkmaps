# Minkowski Maps for a Morphometric Statistical Analysis
#
# Default: NDEBUG
# Debug:   DEBUG=1
#

CXX         = g++
CXXFLAGS    = -Wall -march=native

BOOSTFLAGS  = -Llibdir -L$(BOOST_LIBDIR) -I$(BOOST_INCLUDEDIR) -L$(GSL_LIBDIR) -I$(GSL_INCLUDEDIR)
CXXFLAGS   += $(BOOSTFLAGS)

NDEBUGFLAG  = -O2 -DNDEBUG
DEBUGFLAG   = -D_DEBUG -DDEBUG -O2 -g -pg
ifeq ($(DEBUG),1)
	CXXFLAGS += $(DEBUGFLAG)
else
	CXXFLAGS += $(NDEBUGFLAG)
endif

LDFLAGS     = -lgsl -lgslcblas -lboost_program_options

MINKMAP_GIT_VERSION = $(shell sh -c 'git describe --abbrev=10 --dirty --always')
MINKMAP_WL_OPTIONS = -D__MINKMAP_VERSION=\"$(MINKMAP_GIT_VERSION)\"

# Header
HDR = *.h

# Objects
OBJ = \
	initialize.o \
	helpinghand.o \
	gaussian.o \
	binfield.o \
	randomnumbers.o \
	lookuptable.o \
	lookup_phase_space.o \
	minkowski.o \
	macrostate.o \
	stateclassifier.o \
	binnedskymap.o \
	phasespace.o \

# Binaries
BINARIES = \
	minkmap \

# Default-Target
%.o: %.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $(MOMA_WL_OPTIONS) -o $@ -c $<

minkmap: minkmap.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(MOMA_WL_OPTIONS) -o $@ $(OBJ) $< $(LDFLAGS)

prepare: prepare.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(MOMA_WL_OPTIONS) -o $@ $(OBJ) $< $(LDFLAGS)

all: $(BINARIES)

.PHONY: clean
clean:
	@ rm -rf $(BINARIES) $(OBJ) *~
	@ echo "rm binaries and objects"

