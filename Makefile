# -*- Makefile -*-

### Cori
CXX = CC
CXXFLAGS = -O3 -DMKL_ILP64 -mkl -Wall -Bdynamic -qopenmp
### Sherlock
# CXX = mpiicpc
# CXXFLAGS = -O3 -std=c++14 -DMKL_ILP64 -mkl -Wall -Bdynamic -qopenmp

# MKLLINK = -lpthread -lm -ldl
MKLLINK = -lpthread -lm -ldl
#CXXFLAGS = -O3 -no-prec-div -fp-model fast=2 -xHost -qopenmp
#CXXFLAGS = -O3 -Wall -Bdynamic -fopenmp
#PARPACKLIB = -L/global/homes/t/tatang/lib64
ARPACKHD = -I$(ARPACKINC)
LIBRARY = -L$(ARPACKLIB) 

all:genBasis.out kGroundState.out Spectra.out
.PHONY: all

test: timeTest.out
.PHONY: test

#ED:main.o globalClass.o utils.o
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) -lparpack main.o globalClass.o utils.o -o ED

genBasis.out:genBasisMain.o Geometry.o Basis.o utils.o
	$(CXX) $(CXXFLAGS) $(LIBRARY) genBasisMain.o Geometry.o Basis.o utils.o -o genBasis.out
    
kGroundState.out:kGroundStateMain.o Geometry.o Basis.o Operators.o utils.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) -lparpack kGroundStateMain.o Geometry.o Basis.o Operators.o algebra.o utils.o -o kGroundState.out

Spectra.out:SpectraMain.o Geometry.o Basis.o Operators.o utils.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) -lparpack SpectraMain.o Geometry.o Basis.o Operators.o algebra.o utils.o -o Spectra.out

timeTest.out:timeTest.o algebra.o utils.o Basis.o Geometry.o Operators.o
	$(CXX) $(CXXFLAGS) $(ARPACKHD) $(LIBRARY) -lparpack timeTest.o algebra.o utils.o Basis.o Geometry.o Operators.o -o timeTest.out $(MKLLINK)
#main.o:main.cpp
#	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c main.cpp

genBasisMain.o:genBasisMain.cpp
	$(CXX) $(CXXFLAGS) -c genBasisMain.cpp

kGroundStateMain.o:kGroundStateMain.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c kGroundStateMain.cpp
	
SpectraMain.o:SpectraMain.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c SpectraMain.cpp

timeTest.o:timeTest.cpp
	$(CXX) $(CXXFLAGS) $(ARPACKHD) -c timeTest.cpp

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
