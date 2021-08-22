cmake_minimum_required(VERSION 3.10)

project(ED VERSION 2.0.0 LANGUAGES CXX)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "build on Linux ...")
    set(CMAKE_CXX_COMPILER g++)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(STATUS "build on macOS ...")
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
# gsl
set(GSLDIR $ENV{GSL_DIR})
set(GSL gsl gslcblas)

# Compiler Flags
set(CMAKE_CXX_FLAGS "-O3 ${MKL_FLAGS} -Wall ${OPENMP_FLAG}")
# Include/Library Path
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${ARPACKINC} ${MKLINC} ${BOOSTROOT} ${MPIDIR}/include ${GSLDIR}/include)
link_directories(${ARPACKLIB} ${BOOSTLIB} ${MPIDIR}/lib ${GSLDIR}/lib /opt/intel/lib)
set(LINKS ${MKL} ${BOOST} parpack ${MPI} ${GSL})

add_library(ed STATIC geometry/geometry.cpp basis/basisState.cpp operator/links.cpp measure/config.cpp algebra/algebra.cpp pulse/pulse.cpp utils/comb.cpp utils/paras.cpp utils/path.cpp)
target_link_libraries(ed ${LINKS})

# Excutables
set(EXE_SRC progSrc)
add_executable(main.out ${EXE_SRC}/main.cpp)
target_link_libraries(main.out ed ${LINKS})

# Test
set(TEST_SRC test)
add_executable(test.out ${TEST_SRC}/test.cpp)
target_link_libraries(test.out ed ${LINKS})