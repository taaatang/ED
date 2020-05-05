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
include_directories(Global Geometry Basis Operator Algebra Utils Solver)
include_directories(${ARPACKINC} ${MKLINC})
set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall -qopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${MKL_LINK_FLAGS}")

# SRC and Hearder files
set(GLOBAL_PARA "Global/globalPara.hpp")
set(UTILS "Utils/utils.cpp" "Utils/utils.hpp" "Utils/HelperClass.hpp")
set(GEOMETRY "Geometry/Geometry.cpp" "Geometry/Geometry.hpp")
set(BASIS "Basis/Basis.cpp" "Basis/Basis.hpp")
set(ALGEBRA "Algebra/algebra.cpp" "Algebra/algebra.hpp")
set(OPERATOR "Operator/Operators.cpp" "Operator/Operators.hpp" "Operator/SparseMatrix.hpp")

set(ARPACK "Solver/PARPACKSolver.hpp")
set(LANCZOS "Solver/LANCZOSIterator.hpp")
set(SPECTRA "Solver/Spectra.hpp" ${LANCZOS})
set(TIME "Solver/TimeEvolver.hpp" ${LANCZOS})

set(BASIS_DEP ${BASIS} ${GEOMETRY} ${GLOBAL_PARA} ${UTILS})
set(GROUND_STATE_DEP ${BASIS_DEP} ${OPERATOR} ${ARPACK} ${ALGEBRA})
set(SPECTRA_DEP ${GROUND_STATE_DEP} ${SPECTRA})
set(TIME_DEP ${GROUND_STATE_DEP} ${TIME}) 

# Excutables
set(EXE_SRC ExcutableSRC)

add_executable(genBasis.out ${EXE_SRC}/genBasisMain.cpp ${BASIS_DEP})

add_executable(GroundState.out ${EXE_SRC}/GroundStateMain.cpp ${GROUND_STATE_DEP})
target_link_libraries(GroundState.out ${MKL})
target_link_libraries(GroundState.out PARPACK)

add_executable(Spectra.out ${EXE_SRC}/SpectraMain.cpp ${SPECTRA_DEP})
target_link_libraries(Spectra.out ${MKL})
target_link_libraries(Spectra.out PARPACK)

add_executable(HubbardSpec.out ${EXE_SRC}/HubbardSpecMain.cpp ${SPECTRA_DEP})
target_link_libraries(HubbardSpec.out ${MKL})
target_link_libraries(HubbardSpec.out PARPACK)
                           

# add_executable(Time.out timeTest.cpp algebra.cpp utils.cpp Basis.cpp Geometry.cpp Operators.cpp)
    