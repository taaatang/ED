cmake_minimum_required(VERSION 3.10)

project(ED VERSION 1.0.0 LANGUAGES CXX)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER CC)
set(ARPACKINC $ENV{ARPACKINC})
set(ARPACKLIB $ENV{ARPACKLIB})
set(MKLROOT $ENV{MKLROOT})

set(MKL_FLAGS "-DMKL_ILP64 -I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group")
set(MKL_LINK_FLAGS "-liomp5 -lpthread -lm -ldl")

link_directories(${ARPACKLIB})
include_directories(${ARPACKINC})

set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall -qopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${MKL_LINK_FLAGS}")

# set(GLOBAL "globalPara.hpp")
# set(UTILS "utils.cpp utils.hpp HelperClass.hpp")
# set(GEOMETRY "Geometry.cpp Geometry.hpp")
# set(BASIS "Basis.cpp Basis.hpp")
# set(ALGEBRA "algebra.cpp algebra.hpp")
# set(OPERATOR "Operators.cpp Operators.hpp")

# set(SPARSE "SparseMatrix.hpp")
# set(ARPACK "PARPACKSolver.hpp")
# set(LANCZOS "LANCZOSIterator.hpp")
# set(SPECTRA "Spectra.hpp")
# set(TIME "TimeEvolver.hpp")

# set(BASIS_DEP "${BASIS} ${GEOMETRY} ${GLOBAL} ${UTILS}")
# set(GROUND_STATE_DEP "${BASIS_DEP} ${OPERATOR} ${SPARSE} ${ARPACK} ${ALGEBRA}")
# set(SPECTRA_DEP "${GROUND_STATE_DEP} ${LANCZOS}")
# set(TIME_DEP "${GROUND_STATE_DEP} ${LANCZOS}") 


add_executable(genBasis.out genBasisMain.cpp Geometry.cpp Basis.cpp utils.cpp)
add_executable(kGroundState.out kGroundStateMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
add_executable(Spectra.out SpectraMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
add_executable(HubbardSpec.out HubbardSpecMain.cpp Geometry.cpp Basis.cpp Operators.cpp utils.cpp algebra.cpp)
# add_executable(Time.out timeTest.cpp algebra.cpp utils.cpp Basis.cpp Geometry.cpp Operators.cpp)
    