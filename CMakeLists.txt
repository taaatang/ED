cmake_minimum_required(VERSION 3.10)

set(HOST $ENV{MY_HOST})
if(HOST STREQUAL "cori")
    set(CMAKE_CXX_COMPILER CC)
elseif(HOST STREQUAL "sherlock")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

project(ED VERSION 1.0.0 LANGUAGES CXX)

# PARPACK
set(ARPACKINC $ENV{ARPACKINC})
set(ARPACKLIB $ENV{ARPACKLIB})
set(ARPACK_LINK_FLAGS "-lparpack")
add_library(PARPACK SHARED IMPORTED)
set_target_properties(PARPACK PROPERTIES IMPORTED_LOCATION ${ARPACKLIB}/libparpack.so)

# MKL library with ilp64 support
set(MKLROOT $ENV{MKLROOT})
set(MKLINC ${MKLROOT}/include)
set(MKL_FLAGS "-DMKL_ILP64")
set(MKL_LINK_FLAGS "-liomp5 -lpthread -lm -ldl")
add_library(MKL_ilp64 STATIC IMPORTED)
set_target_properties(MKL_ilp64 PROPERTIES IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a)
add_library(MKL_thread STATIC IMPORTED)
set_target_properties(MKL_thread PROPERTIES IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_intel_thread.a)
add_library(MKL_core STATIC IMPORTED)
set_target_properties(MKL_core PROPERTIES IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_core.a) 
set(MKL -Wl,--start-group MKL_ilp64 MKL_thread MKL_core -Wl,--end-group)

# Include/Library Path, Compiler/Linker Flags
include_directories(${ARPACKINC} ${MKLINC})
set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall -qopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${MKL_LINK_FLAGS}")

# SRC and Hearder files
set(GLOBAL "globalPara.hpp")
set(UTILS "utils.cpp" "utils.hpp" "HelperClass.hpp")
set(GEOMETRY "Geometry.cpp" "Geometry.hpp")
set(BASIS "Basis.cpp" "Basis.hpp")
set(ALGEBRA "algebra.cpp" "algebra.hpp")
set(OPERATOR "Operators.cpp" "Operators.hpp")

set(SPARSE "SparseMatrix.hpp")
set(ARPACK "PARPACKSolver.hpp")
set(LANCZOS "LANCZOSIterator.hpp")
set(SPECTRA "Spectra.hpp")
set(TIME "TimeEvolver.hpp")

set(BASIS_DEP ${BASIS} ${GEOMETRY} ${GLOBAL} ${UTILS})
set(GROUND_STATE_DEP ${BASIS_DEP} ${OPERATOR} ${SPARSE} ${ARPACK} ${ALGEBRA})
set(SPECTRA_DEP ${GROUND_STATE_DEP} ${LANCZOS})
set(TIME_DEP ${GROUND_STATE_DEP} ${LANCZOS} ${TIME}) 

# Excutables
add_executable(genBasis.out genBasisMain.cpp ${BASIS_DEP})

add_executable(kGroundState.out kGroundStateMain.cpp ${GROUND_STATE_DEP})
target_link_libraries(kGroundState.out ${MKL})
target_link_libraries(kGroundState.out PARPACK)

add_executable(Spectra.out SpectraMain.cpp ${SPECTRA_DEP})
target_link_libraries(Spectra.out ${MKL})
target_link_libraries(Spectra.out PARPACK)

add_executable(HubbardSpec.out HubbardSpecMain.cpp ${SPECTRA_DEP})
target_link_libraries(HubbardSpec.out ${MKL})
target_link_libraries(HubbardSpec.out PARPACK)
                           

# add_executable(Time.out timeTest.cpp algebra.cpp utils.cpp Basis.cpp Geometry.cpp Operators.cpp)
    