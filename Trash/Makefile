# -*- Makefile -*-

MKLFLAG = -DMKL_ILP64 -I${MKLROOT}/include
# --start-group ... --end-group: resolving circular dependences between several libraries. -Wl:pass options to the linker
MKLLINK = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

### Cori
CXX = CC
CXXFLAGS = -O3 $(MKLFLAG) -w2 -qopenmp

### Sherlock
# CXX = mpiicpc
# CXXFLAGS = -O3 -std=c++14 $(MKLFLAG) -Wall -Bdynamic -qopenmp

# CXXFLAGS = -O3 -std=c++14 -DMKL_ILP64 -mkl -Wall -Bdynamic -qopenmp
#CXXFLAGS = -O3 -no-prec-div -fp-model fast=2 -xHost -qopenmp
#CXXFLAGS = -O3 -Wall -Bdynamic -fopenmp

ARPACKHD = -I$(ARPACKINC)
LIBRARY = -L$(ARPACKLIB)

all:genBasis.out kGroundState.out Spectra.out HubbardSpec.out
.PHONY: all

test: timeTest.out geometryTest.out
.PHONY: test

#ED:main.o globalClass.o utils.o
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) -lparpack main.o globalClass.o utils.o -o ED

genBasis.out:genBasisMain.o Geometry.o Basis.o utils.o
	$(CXX) $(CXXFLAGS) $(LIBRARY) genBasisMain.o Geometry.o Basis.o utils.o -o genBasis.out
    
kGroundState.out:kGroundStateMain.o Geometry.o Basis.o Operators.o utils.o algebra.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY)  kGroundStateMain.o Geometry.o Basis.o Operators.o algebra.o utils.o -o kGroundState.out $(MKLLINK) -lparpack

Spectra.out:SpectraMain.o Geometry.o Basis.o Operators.o utils.o algebra.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY)  SpectraMain.o Geometry.o Basis.o Operators.o algebra.o utils.o -o Spectra.out $(MKLLINK) -lparpack
HubbardSpec.out:HubbardSpecMain.o Geometry.o Basis.o Operators.o utils.o algebra.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY)  HubbardSpecMain.o Geometry.o Basis.o Operators.o algebra.o utils.o -o HubbardSpec.out $(MKLLINK) -lparpack

timeTest.out:timeTest.o algebra.o utils.o Basis.o Geometry.o Operators.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) timeTest.o algebra.o utils.o Basis.o Geometry.o Operators.o -o timeTest.out $(MKLLINK) -lparpack
geometryTest.out:geometryTest.o Geometry.o utils.o
	$(CXX) $(CXXFLAGS) geometryTest.o Geometry.o utils.o -o geometryTest.out
#main.o:main.cpp
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c main.cpp

genBasisMain.o:genBasisMain.cpp
	$(CXX) $(CXXFLAGS) -c genBasisMain.cpp

kGroundStateMain.o:kGroundStateMain.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c kGroundStateMain.cpp
	
SpectraMain.o:SpectraMain.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c SpectraMain.cpp

HubbardSpecMain.o:HubbardSpecMain.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c HubbardSpecMain.cpp

timeTest.o:timeTest.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c timeTest.cpp
geometryTest.o:geometryTest.cpp
	$(CXX) $(CXXFLAGS) -c geometryTest.cpp

#globalClass.o:globalClass.cpp
#	$(CXX) $(CXXFLAGS) -c globalClass.cpp

Geometry.o:Geometry.cpp
	$(CXX) $(CXXFLAGS) -c Geometry.cpp

Basis.o:Basis.cpp
	$(CXX) $(CXXFLAGS) -c Basis.cpp

Operators.o:Operators.cpp
	$(CXX) $(CXXFLAGS) -c Operators.cpp

utils.o:utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

algebra.o:algebra.cpp
	$(CXX) $(CXXFLAGS) -c algebra.cpp
    
#TEST:test.o
#	$(CXX) $(CXXFLAGS) $(LIBRARY) -larpack test.o -o TEST
#test.o: test.cpp
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c test.cpp
#
#PTEST:ptest.o
#	$(CXX) $(CXXFLAGS) $(LIBRARY) -lparpack ptest.o -o PTEST
#ptest.o: ptest.cpp
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c ptest.cpp

clean:
	@echo "clean up..."
	rm -f *.o
	rm -f *.out
.PHONY: clean
