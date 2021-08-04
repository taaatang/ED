cmake_minimum_required(VERSION 3.10)

set(HOST $ENV{MY_HOST})
if(HOST STREQUAL "cori")
    set(CMAKE_CXX_COMPILER CC)
    set(OPENMP_FLAG -qopenmp)
elseif(HOST STREQUAL "sherlock" OR HOST STREQUAL "acer")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(OPENMP_FLAG -qopenmp)
elseif(HOST STREQUAL "mac")
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-11)
    set(OPENMP_FLAG -fopenmp)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(ED VERSION 1.0.0 LANGUAGES CXX)

# PARPACK
set(ARPACKINC $ENV{ARPACKINC})
set(ARPACKLIB $ENV{ARPACKLIB})

# MKL library with ilp64 support
set(MKLROOT $ENV{MKLROOT})
set(MKLINC ${MKLROOT}/include)

if(HOST STREQUAL "mac")
    set(MPIINC /usr/local/Cellar/open-mpi/4.1.1_2/include)
    set(MPILIB /usr/local/Cellar/open-mpi/4.1.1_2/lib)
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
if(HOST STREQUAL "cori")
    include_directories(/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/include)
    link_directories(/global/common/sw/cray/cnl7/haswell/gsl/2.5/intel/19.0.3.199/7twqxxq/lib)
endif()


# SRC and Hearder files
set(GLOBAL_PARA "Global/globalType.hpp" "Global/globalPara.hpp" "Global/config.hpp")
set(UTILS "Utils/comb.cpp" "Utils/io.cpp"  "Utils/paras.cpp" "Utils/bitop.hpp" "Utils/comb.hpp" "Utils/io.hpp" "Utils/mpiwrap.hpp" "Utils/paras.hpp" "Utils/path.hpp" "Utils/timer.hpp")
set(PULSE "Pulse/pulse.cpp" "Pulse/pulse.hpp")
set(GEOMETRY "Geometry/Geometry.cpp" "Geometry/Geometry.hpp" "Geometry/Vec3.hpp" "Geometry/Transform.hpp" "Geometry/Generator.hpp")
set(BASIS "Basis/BasisState.cpp" "Basis/BasisState.hpp" "Basis/Basis.hpp")
set(ALGEBRA "Algebra/algebra.cpp" "Algebra/algebra.hpp")
set(OPERATOR "Operator/links.cpp" "Operator/Operators.cpp" "Operator/localOperators.hpp" "Operator/Operations.hpp" "Operator/OperatorsBase.hpp"  "Operator/Operators.hpp" "Operator/SparseMatrix.hpp" "Operator/links.hpp")

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
add_executable(test.out ${TEST_SRC}/basisState.cpp   ${BASIS} ${GEOMETRY} ${GLOBAL_PARA} Operator/localOperators.hpp "Utils/comb.cpp" "Utils/io.cpp"  "Utils/bitop.hpp" "Utils/comb.hpp" "Utils/io.hpp" "Utils/mpiwrap.hpp" )
target_link_libraries(test.out ${LINKS} gsl gslcblas)

