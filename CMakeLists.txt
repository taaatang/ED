cmake_minimum_required(VERSION 3.10)

project(ED VERSION 1.0.0 LANGUAGE CXX)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

set(MKL_FLAGS "-DMKL_ILP64 -I$ENV(MKLROOT)/include")
set(MKL_LINK_FLAGS "-Wl,--start-group $ENV(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $ENV(MKLROOT)/lib/intel64/libmkl_intel_thread.a $ENV(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl")

set(ARPACKINC $ENV(ARPACKINC))
set(ARPACKLIB $ENV(ARPACKLIB))

link_directories(${ARPACKLIB})
include_directories(${ARPACKINC})

set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -W2 -qopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${MKL_LINK_FLAGS}")

add_executable(genBasis.out genBasisMain.cpp Geometry.cpp Basis.o utils.cpp)
add_executable(kGroundState.out kGroundStateMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
add_executable(Spectra.out SpectraMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
add_executable(HubbardSpec.out HubbardSpecMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
# add_executable(Time.out timeTest.cpp algebra.cpp utils.cpp Basis.cpp Geometry.cpp Operators.cpp)
    