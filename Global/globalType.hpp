#ifndef __GLOBALTYPE_H__
#define __GLOBALTYPE_H__

#include <iostream>
#include <cmath>
#include <assert.h>
#include <chrono>
#include <vector>
#include <map>
#include <unordered_map>
#include <complex>
#include <string>
#include <utility>
#include <omp.h>

// for filesystem. if not c++ 17, use boost library
#define CPP_17
/*
    *****************
    * Compile Macro *
    * ***************
*/

// #define TEST
// #undef TEST

#define OMP_
// #undef OMP_

#define SAXPY
#undef SAXPY

/*
    **************************
    * Hamiltonian Model Type *
    **************************
*/
/*
    D4m: for multi-band hubbard, due to orbital phases, the symmetry is only a sub group {I, RZ}. R is rotation, Z is reflection.
    D4m5: include 2 Pz orbitals. the symmetry is {I, RZ}x{I, M}. M is the mirror symmetry about the xy plane.
*/
enum PointGroup {NONE,D6, C6, D4, D4m, D4m5, C4, D3, C3};
enum ORBITAL {SINGLE, Dx2y2, Dz2, Px, Py, Pzu, Pzd, TypeA, TypeB};
enum LATTICE_MODEL {HUBBARD,tJ,HEISENBERG};
enum LINK_TYPE {SUPER_EXCHANGE_J, CHIRAL_K, HOPPING_T, CHARGE_TRANSFER_V, HUBBARD_U, EXCHANGE_J, PAIR_HOPPING_J};

/*
    *************************
    * Hamiltonian Data Type *
    *************************
*/

#define COMPLEX_OP // complex valued operator (Hq, Szq ...)
#define COMPLEX_DATA // set dataType = cdouble
// #undef COMPLEX_DATA // det dataType = double

/*
    ************
    * Symmetry *
    ************
*/

#define TRANSLATION_SYMM

/*
    ***************
    * Type Define *
    ***************
*/
// Matrix index datatype
typedef unsigned long long int idx_t; 
typedef unsigned int uint;
typedef std::complex<double> cdouble;
// #ifdef COMPLEX_DATA
typedef cdouble dataType;
// #else
//     typedef double dataType;
// #endif

typedef std::unordered_map<idx_t,dataType> MAP;
typedef std::unordered_map<idx_t,dataType>::iterator MAPIT;

/*
    ******************
    * Geometry Class *
    ******************
*/
typedef int BasisLat[2];
typedef double BasisXY[2];
typedef std::vector<double> VecD;
typedef std::vector<int> VecI;
typedef std::vector<std::vector<int>> BondMap;

/*
    ***************
    * Basis Class *
    ***************
*/
// use the binary form of an integer to represent basis state
#define BINARY_REP
// basis distributed among MPI workers
// #define DISTRIBUTED_BASIS

#ifdef DISTRIBUTED_BASIS
    // wavefunction vectors distributed among MPI workers
    #define DISTRIBUTED_STATE
#endif

typedef std::pair<idx_t,idx_t> pairIndex;
#define KEEP_BASIS_NORM

/*
    ******************
    * Spin Operators *
    ******************
*/

enum SPIN {SPIN_UP, SPIN_DOWN, SPIN_UD};


/*
    **********************************
    * Sparse Matrix PARTITION Method *
    **********************************
*/

enum MATRIX_PARTITION {ROW_PARTITION, COL_PARTITION};

/*
    *******************
    * LANCZOSIterator *
    *******************
*/

enum LANCZOS_OPTION {ALPHA_BETA, ALPHA_BETA_Q};

#endif // __GLOBALTYPE_H__