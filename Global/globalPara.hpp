//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
#ifndef globalPara_hpp
#define globalPara_hpp

#include <iostream>
#include <assert.h>
#include <chrono>
#include <vector>
#include <map>
#include <unordered_map>
#include <complex>
#include <string>
#include <utility>
#include <omp.h>

/*
    *****************
    * Compile Macro *
    * ***************
*/

// #define TEST
// #undef TEST

#define OMP_
// #undef OMP_

const bool DEBUG = false;
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
enum ORBITAL {SINGLE, Dx2y2, Px, Py, Pzu, Pzd, TypeA, TypeB};
enum LATTICE_MODEL {HUBBARD,t_J,HEISENBERG};
enum LINK_TYPE {SUPER_EXCHANGE_J, CHIRAL_K, HOPPING_T, CHARGE_TRANSFER_V, HUBBARD_U};
const LATTICE_MODEL MODEL=HEISENBERG;

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
    *******
    * MPI *
    *******
*/
const int MPI_MASTER = 0;

/*
    *****************
    * File I/O Path *
    *****************
*/
// cori
// const std::string ROOT_DATA_PATH = "/global/project/projectdirs/m2757/tatang";  
// sherlock
// const std::string ROOT_DATA_PATH = "/oak/stanford/groups/tpd/tatang"; 
// acer
const std::string ROOT_DATA_PATH = "/home/tatang/Project/DATA"; 

// const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/Photodoping";
// const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/TriAngtJ";
const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/TriAngHeis";
/*
    ******************
    * Math Constants *
    ******************
*/

const double INFINITESIMAL = 1e-8;
const double PI = 3.14159265358979323846;
const std::complex<double> CPLX_I(0.0, 1.0);

/*
    ***************
    * Type Define *
    ***************
*/
// Matrix index datatype
typedef unsigned long long int ind_int; 
typedef unsigned int uint;
typedef std::complex<double> cdouble;
// #ifdef COMPLEX_DATA
typedef std::complex<double> dataType;
// #else
//     typedef double dataType;
// #endif

typedef std::unordered_map<ind_int,dataType> MAP;
typedef std::unordered_map<ind_int,dataType>::iterator MAPIT;

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
const int MAX_BOND_NUM = 6;

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

typedef std::pair<ind_int,ind_int> pairIndex;
#define KEEP_BASIS_NORM

/*
    ******************
    * Spin Operators *
    ******************
*/

const int SPIN_DIM = 2;
enum SPIN {SPIN_UP, SPIN_DOWN, SPIN_UD};


/*
    **********************************
    * Sparse Matrix PARTITION Method *
    **********************************
*/

enum MATRIX_PARTITION {ROW_PARTITION, COL_PARTITION};
#ifdef DISTRIBUTED_STATE
    #define SPM_COL_PARTITION
    const MATRIX_PARTITION PARTITION = COL_PARTITION;
#else
    #define SPM_ROW_PARTITION
    const MATRIX_PARTITION PARTITION = ROW_PARTITION;
#endif


/*
    *****************
    * PARPACKSolver *
    *****************
*/
// MAX parpack iterations
const int PARPACK_MAXITERATION = 5000;
//MIN and MAX ncv
const int PARPACK_MINNCV = 7;
// const int PARPACK_MAXNCV = 500;

/*
    *******************
    * LANCZOSIterator *
    *******************
*/

enum LANCZOS_OPTION {ALPHA_BETA, ALPHA_BETA_Q};
const LANCZOS_OPTION LANCZOS_DEFAULT_OPTION = ALPHA_BETA;

#endif // globalPara_hpp
