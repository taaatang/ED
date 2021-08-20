cmake_minimum_required(VERSION 3.10)

project(ED VERSION 1.0.0 LANGUAGES CXX)

# set(HOST $ENV{MY_HOST})
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_COMPILER g++)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-11)
endif()
set(OPENMP_FLAG -fopenmp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# PARPACK
set(ARPACKINC $ENV{ARPACKINC})
set(ARPACKLIB $ENV{ARPACKLIB})

# MKL library with ilp64 support
set(MKLROOT $ENV{MKLROOT})
set(MKLINC ${MKLROOT}/include)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(MKLLIB ${MKLROOT}/lib/intel64)
    set(MPIDIR $ENV{MPICH_DIR})
    set(MPI mpich)
else(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(MKLLIB ${MKLROOT}/lib)
    set(MPIDIR $ENV{OMPI_DIR})
    set(MPI mpi)
endif()

set(MKL_FLAGS "-DMKL_ILP64")
add_library(MKL_ilp64 STATIC IMPORTED)
set_target_properties(MKL_ilp64 PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_intel_ilp64.a)

add_library(MKL_thread STATIC IMPORTED)
set_target_properties(MKL_thread PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_intel_thread.a)

add_library(MKL_core STATIC IMPORTED)
set_target_properties(MKL_core PROPERTIES IMPORTED_LOCATION ${MKLLIB}/libmkl_core.a)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(MKL -Wl,--start-group MKL_ilp64 MKL_thread MKL_core -Wl,--end-group iomp5 pthread m dl)
else()
    set(MKL MKL_ilp64 MKL_thread MKL_core MKL_ilp64 MKL_thread MKL_core iomp5 pthread m dl)
endif()
# Boost
set(BOOSTROOT $ENV{BOOST_ROOT})
set(BOOSTLIB $ENV{BOOST_LIB})
set(BOOST boost_system boost_filesystem)

set(GSLDIR $ENV{GSL_DIR})
set(GSL gsl gslcblas)

# Compiler Flags
set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall ${OPENMP_FLAG}")
# Include/Library Path
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${ARPACKINC} ${MKLINC} ${BOOSTROOT} ${MPIDIR}/include ${GSLDIR}/include)
link_directories(${ARPACKLIB} ${BOOSTLIB} ${MPIDIR}/lib ${GSLDIR}/lib /opt/intel/lib)
set(LINKS ${MKL} ${BOOST} parpack ${MPI} ${GSL})

# SRC and Hearder files
set(GLOBAL_PARA "Global/globalType.hpp" "Global/globalPara.hpp")
set(UTILS "Utils/comb.cpp" "Utils/io.cpp"  "Utils/paras.cpp" "Utils/bitop.hpp" "Utils/comb.hpp" "Utils/io.hpp" "Utils/mpiwrap.hpp" "Utils/paras.hpp" "Utils/path.hpp" "Utils/timer.hpp")
set(PULSE "Pulse/pulse.cpp" "Pulse/pulse.hpp")
set(GEOMETRY "Geometry/Geometry.cpp" "Geometry/Geometry.hpp" "Geometry/Vec3.hpp" "Geometry/Transform.hpp" "Geometry/Generator.hpp")
set(BASIS "Basis/BasisState.cpp" "Basis/BasisState.hpp" "Basis/Basis.hpp")
set(ALGEBRA "Algebra/algebra.cpp" "Algebra/algebra.hpp")
set(OPERATOR "Operator/links.cpp" "Operator/localOperators.hpp" "Operator/Operations.hpp" "Operator/OperatorsBase.hpp"  "Operator/Operators.hpp" "Operator/SparseMatrix.hpp" "Operator/links.hpp")

set(ARPACK "Solver/PARPACKSolver.hpp")
set(LANCZOS "Solver/LANCZOSIterator.hpp")
set(BiCGSTAB "Solver/BiCGSTAB.hpp")
set(SPECTRA "Solver/Spectra.hpp" ${LANCZOS} ${BiCGSTAB})
set(TIME "Solver/TimeEvolver.hpp" ${LANCZOS})

set(GEOMETRY_DEP ${GEOMETRY} ${GLOBAL_PARA} ${UTILS})
set(BASIS_DEP ${BASIS} ${GEOMETRY_DEP})
set(GROUND_STATE_DEP ${BASIS_DEP} ${OPERATOR} ${ARPACK} ${ALGEBRA})
set(SPECTRA_DEP ${GROUND_STATE_DEP} ${SPECTRA})
set(TIME_DEP ${GROUND_STATE_DEP} ${TIME} ${PULSE})
set(ALL_DEP ${SPECTRA_DEP} ${TIME} ${PULSE})

# Excutables
set(EXE_SRC ExcutableSRC)

add_executable(main.out ${EXE_SRC}/main.cpp ${ALL_DEP})
target_link_libraries(main.out ${LINKS} gsl gslcblas)

# Test
set(TEST_SRC Test)

#add_executable(test.out ${TEST_SRC}/basisState.cpp ${ALL_DEP})
add_executable(test.out ${TEST_SRC}/test.cpp ${ALGEBRA} ${PULSE} ${BASIS} ${GEOMETRY} ${GLOBAL_PARA} Operator/localOperators.hpp "Utils/comb.cpp" "Utils/io.cpp"  "Utils/bitop.hpp" "Utils/comb.hpp" "Utils/io.hpp" "Utils/mpiwrap.hpp" "Utils/paras.cpp" "Utils/paras.hpp" "Utils/path.cpp" "Utils/path.hpp" "Pulse/pulse.cpp" "Pulse/pulse.hpp" "Measure/config.cpp" "Measure/config.hpp")
target_link_libraries(test.out ${LINKS} gsl gslcblas)

