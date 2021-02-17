cmake_minimum_required(VERSION 3.10)

set(HOST $ENV{MY_HOST})
if(HOST STREQUAL "cori")
    set(CMAKE_CXX_COMPILER CC)
    set(OPENMP_FLAG -qopenmp)
elseif(HOST STREQUAL "sherlock" OR HOST STREQUAL "acer")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
    set(OPENMP_FLAG -qopenmp)
elseif(HOST STREQUAL "mac")
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-10)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
    set(OPENMP_FLAG -fopenmp)
endif()

project(ED VERSION 1.0.0 LANGUAGES CXX)

# PARPACK
set(ARPACKINC $ENV{ARPACKINC})
set(ARPACKLIB $ENV{ARPACKLIB})

# MKL library with ilp64 support
set(MKLROOT $ENV{MKLROOT})
set(MKLINC ${MKLROOT}/include)

if(HOST STREQUAL "mac")
    set(MPIINC /usr/local/Cellar/open-mpi/4.0.5/include)
    set(MPILIB /usr/local/Cellar/open-mpi/4.0.5/lib)
    set(MPI mpi)
    set(MKLLIB ${MKLROOT}/lib)
else()
    set(MPIINC )
    set(MPILIB ) 
    set(MPI )
    set(MKLLIB ${MKLROOT}/lib/intel64)
endif()

set(MKL_FLAGS "-DMKL_ILP64")
add_library(MKL_ilp64 STATIC IMPORTED)
set_target_properties(MKL_ilp64 PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_intel_ilp64.a)
add_library(MKL_thread STATIC IMPORTED)
set_target_properties(MKL_thread PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_intel_thread.a)
add_library(MKL_core STATIC IMPORTED)
set_target_properties(MKL_core PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_core.a) 
if(HOST STREQUAL "mac")
    set(MKL MKL_ilp64 MKL_thread MKL_core MKL_ilp64 MKL_thread MKL_core iomp5 pthread m dl)
else()
    set(MKL -Wl,--start-group MKL_ilp64 MKL_thread MKL_core -Wl,--end-group iomp5 pthread m dl)
endif()
# Boost
set(BOOSTROOT $ENV{BOOST_ROOT})
set(BOOSTLIB $ENV{BOOST_LIB})
set(BOOST boost_system boost_filesystem)

set(LINKS ${MKL} ${BOOST} parpack ${MPI})

# Compiler Flags
set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall ${OPENMP_FLAG}")
# Include/Library Path
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${ARPACKINC} ${MKLINC} ${BOOSTROOT} ${MPIINC})
link_directories(${ARPACKLIB} ${BOOSTLIB} ${MPILIB} /opt/intel/lib)
if(HOST STREQUAL "mac")
    include_directories(/usr/local/Cellar/gsl/2.6/include)
    link_directories(/usr/local/Cellar/gsl/2.6/lib)
endif()


# SRC and Hearder files
set(GLOBAL_PARA "Global/globalType.hpp" "Global/globalPara.hpp" "Global/config.hpp")
set(UTILS "Utils/utils.cpp" "Utils/paras.cpp" "Utils/utils.hpp" "Utils/timer.hpp" "Utils/paras.hpp")
set(PULSE "Pulse/pulse.cpp" "Pulse/pulse.hpp")
set(GEOMETRY "Geometry/Geometry.cpp" "Geometry/Geometry.hpp")
set(BASIS "Basis/Basis.cpp" "Basis/Basis.hpp")
set(ALGEBRA "Algebra/algebra.cpp" "Algebra/algebra.hpp")
set(OPERATOR "Operator/links.cpp" "Operator/Operators.cpp" "Operator/Operations.hpp" "Operator/OperatorsBase.hpp"  "Operator/Operators.hpp" "Operator/SparseMatrix.hpp" "Operator/links.hpp")

set(ARPACK "Solver/PARPACKSolver.hpp")
set(LANCZOS "Solver/LANCZOSIterator.hpp")
set(SPECTRA "Solver/Spectra.hpp" ${LANCZOS})
set(TIME "Solver/TimeEvolver.hpp" ${LANCZOS})

set(GEOMETRY_DEP ${GEOMETRY} ${GLOBAL_PARA} ${UTILS})
set(BASIS_DEP ${BASIS} ${GEOMETRY_DEP})
set(GROUND_STATE_DEP ${BASIS_DEP} ${OPERATOR} ${ARPACK} ${ALGEBRA})
set(SPECTRA_DEP ${GROUND_STATE_DEP} ${SPECTRA})
set(TIME_DEP ${GROUND_STATE_DEP} ${TIME} ${PULSE})

# Excutables
set(EXE_SRC ExcutableSRC)

add_executable(main.out ${EXE_SRC}/main.cpp ${TIME_DEP})
target_link_libraries(main.out ${LINKS} gsl gslcblas)

add_executable(genBasis.out ${EXE_SRC}/genBasisMain.cpp ${BASIS_DEP})
target_link_libraries(genBasis.out ${BOOST})

add_executable(GroundState.out ${EXE_SRC}/GroundStateMain.cpp ${GROUND_STATE_DEP})
target_link_libraries(GroundState.out ${LINKS})

add_executable(Spectra.out ${EXE_SRC}/SpectraMain.cpp ${SPECTRA_DEP})
target_link_libraries(Spectra.out ${LINKS})

add_executable(HubbardSpec.out ${EXE_SRC}/HubbardSpecMain.cpp ${SPECTRA_DEP})
target_link_libraries(HubbardSpec.out ${LINKS})

add_executable(HubbardCond.out ${EXE_SRC}/HubbardCondMain.cpp ${SPECTRA_DEP})
target_link_libraries(HubbardCond.out ${LINKS})

# Test
set(TEST_SRC Test)

add_executable(geometryTest.out ${TEST_SRC}/geometryTest.cpp ${GEOMETRY_DEP})

add_executable(timeTest.out ${TEST_SRC}/timeTest.cpp ${TIME_DEP})
target_link_libraries(timeTest.out ${LINKS} parpack gsl gslcblas)

# add_executable(Time.out timeTest.cpp algebra.cpp utils.cpp Basis.cpp Geometry.cpp Operators.cpp) 
