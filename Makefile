UNAME:=$(shell uname -as | tr a-z A-Z)
UNAMEARG:=$(shell uname -s | tr a-z A-Z)
PED_DEPTH_DYLIB_NAME:=costFunction

CLANG_PATH=g++
HOSTNAME= $(shell hostname)
ifneq (,$(findstring analyzethis, $(HOSTNAME) ))
	CLANG_PATH=/opt/local/bin/clang++-mp-3.3 
endif

INCLUDE=

ifneq (,$(findstring DARWIN, $(UNAME)))
	INCLUDE = -I.
	BOOST_PATH=/opt/local
	LDFLAGS=$(BOOST_PATH)/lib/libboost_program_options-mt.a
	#CC=g++
	CC=${CLANG_PATH}
	CFLAGS=$(INCLUDE) -D$(UNAMEARG) -DUSE_PTHREADS -DUSE_LIBDISPATCH -static -Weverything -Werror -Wno-padded -pedantic -Wno-c++11-extensions -Wno-c++98-compat 
	MAKE=make
	EXEC_NAME=simworm_darwin
	DYLIB_SUFFIX=.dylib
endif
ifneq (,$(findstring FREEBSD, $(UNAME)))
	INCLUDE = -I.
	BOOST_PATH=/opt/local
	LDFLAGS=-lboost_program_options -ldispatch -L/usr/lib -L/usr/local/lib -fPIC -lc++ -lm -lBlocksRuntime
	CC=clang
	CFLAGS=$(INCLUDE) -D$(UNAMEARG) -fPIC -fblocks -Weverything -Werror -Wno-padded -pedantic -Wno-c++11-extensions -Wno-c++98-compat -DUSE_LIBDISPATCH
        # CFLAGS=$(INCLUDE) -D$(UNAMEARG) -fPIC -fblocks -Weverything -Werror -Wno-padded -pedantic -Wno-c++11-extensions -Wno-c++98-compat 
	MAKE=gmake
	EXEC_NAME=simworm_bsd
	DYLIB_SUFFIX=.so
endif
ifneq (,$(findstring LINUX, $(UNAME)))
	ifneq (,$(findstring BRUTUS2L, $(UNAME)))
		BOOST_PATH=/usr
		LDFLAGS=/usr/lib64/libboost_program_options.so
		CC=clang
		CFLAGS=$(INCLUDE) -D$(UNAMEARG) -static
		MAKE=make
        	EXEC_NAME=simworm_linux_brutus2l
		DYLIB_SUFFIX=_FC16.so
	endif
	ifneq (,$(findstring HPC, $(UNAME)))
		BOOST_PATH=/data/users/mchiang3/boost_1_47_0/
		LDFLAGS=-fPIC /data/users/mchiang3/boost_1_47_0/stage/lib/libboost_program_options.so
		CC=clang
		CFLAGS=$(INCLUDE) -fPIC -D$(UNAMEARG) -DUSE_PTHREADS 
		MAKE=make
        	EXEC_NAME=simworm_linux_oit
		DYLIB_SUFFIX=_linux.so
	endif		
endif

ifneq (,$(findstring analyzethis, $(HOSTNAME) ))
        CFLAGS+= -Wno-sign-conversion
endif


INCLUDE += -isystem $(BOOST_PATH) -isystem $(BOOST_PATH)/include -isystem Eigen
DEBUG_FLAG=-g
RELEASE_FLAG=-g -O3
MODE=

all: release
debug: CFLAGS+=$(DEBUG_FLAG)
debug: MODE=debug
debug: preprocess_debug
debug: $(EXEC_NAME) dylibs
release: CFLAGS+=$(RELEASE_FLAG)
release: MODE=release
release: preprocess_release
release: $(EXEC_NAME) dylibs
clean:
	rm -f *.o $(EXEC_NAME)
	rm -f *.o $(PED_DEPTH_DYLIB_NAME)$(DYLIB_SUFFIX)
	rm -f gitVersion.h
	rm -f simworm_ccfit
preprocess_release:
	@echo Building $(EXEC_NAME), release build for $(UNAMEARG)...
preprocess_debug:
	@echo Building $(EXEC_NAME), debug build for $(UNAMEARG)...


ProgramOptionParser.o: ProgramOptionParser.cpp CellCycleProfile.h GeometryProfile.h
	$(CC) $(CFLAGS) -c ProgramOptionParser.cpp
CellCycleProfile.o: CellCycleProfile.cpp util.h Array1D.h Array2D.h
	$(CC) $(CFLAGS) -c CellCycleProfile.cpp
MeioticFractionProfile.o: MeioticFractionProfile.cpp Array1D.h Array2D.h
	$(CC) $(CFLAGS) -c MeioticFractionProfile.cpp
GeometryProfile.o: GeometryProfile.cpp util.h Array1D.h Array2D.h definitions.h
	$(CC) $(CFLAGS) -c GeometryProfile.cpp
GermCell.o: GermCell.cpp util.h Array2D.h
	$(CC) $(CFLAGS) -c GermCell.cpp
GermLine.o: GermLine.cpp util.h GeometryProfile.h CellCycleProfile.h Array2D.h GermCell.h ProgramOptionParser.h MeioticFractionProfile.h
	$(CC) $(CFLAGS) -c GermLine.cpp
pedigree_depth_optimization.o: pedigree_depth_optimization.cpp definitions.h costFunction.h
	$(CC) $(CFLAGS) -c pedigree_depth_optimization.cpp
costFunction.o: costFunction.cpp definitions.h ProgramOptionParser.h util.h GermLine.h
	$(CC) $(CFLAGS) -c costFunction.cpp
util.o: util.cpp
	$(CC) $(CFLAGS) -c util.cpp
FitData.o: FitData.cpp
	$(CC) $(CFLAGS) -c FitData.cpp
PhaseIndexProfile.o: PhaseIndexProfile.cpp
	$(CC) $(CFLAGS) -c PhaseIndexProfile.cpp
util_germline.o: util_germline.cpp
	$(CC) $(CFLAGS) -c util_germline.cpp
cell_cycle_fit.o: cell_cycle_fit.cpp ccfit.h
	$(CC) $(CFLAGS) -c cell_cycle_fit.cpp
ccfit.o: ccfit.cpp
	$(CC) $(CFLAGS) -c ccfit.cpp

FORCE:
gitVersion.h: FORCE
	echo "static const char *gitVersion = \"$(shell git rev-parse HEAD && git diff-index HEAD)\";" > $@
run.o: run.cpp costFunction.h gitVersion.h
	$(CC) $(CFLAGS) -c run.cpp


$(EXEC_NAME): FitData.o ccfit.o cell_cycle_fit.o PhaseIndexProfile.o util_germline.o ProgramOptionParser.o MeioticFractionProfile.o CellCycleProfile.o GeometryProfile.o GermCell.o GermLine.o util.o costFunction.o pedigree_depth_optimization.o
	$(CC) FitData.o PhaseIndexProfile.o util_germline.o ProgramOptionParser.o MeioticFractionProfile.o CellCycleProfile.o GeometryProfile.o GermCell.o GermLine.o util.o costFunction.o pedigree_depth_optimization.o -o $(EXEC_NAME) $(LDFLAGS)
	$(CC) FitData.o PhaseIndexProfile.o util_germline.o ProgramOptionParser.o MeioticFractionProfile.o CellCycleProfile.o GeometryProfile.o GermCell.o GermLine.o util.o costFunction.o ccfit.o cell_cycle_fit.o -o simworm_ccfit $(LDFLAGS)
dylibs: PhaseIndexProfile.o util_germline.o ProgramOptionParser.o MeioticFractionProfile.o CellCycleProfile.o GeometryProfile.o GermCell.o GermLine.o util.o costFunction.o run.o
	$(CC) -shared run.o FitData.o PhaseIndexProfile.o util_germline.o ProgramOptionParser.o MeioticFractionProfile.o CellCycleProfile.o GeometryProfile.o GermCell.o GermLine.o util.o costFunction.o -o $(PED_DEPTH_DYLIB_NAME)$(DYLIB_SUFFIX) $(LDFLAGS)
