#ifndef __GLOBALTYPE_H__
#define __GLOBALTYPE_H__

#include <vector>
#include <complex>
#include <utility> // pair
#include <unordered_map>

enum class LATTICE_MODEL {HUBBARD,tJ,HEISENBERG};
/*
    D4m: for multi-band hubbard, due to orbital phases, the symmetry is only a sub group {I, RZ}. R is rotation, Z is reflection.
    D4m5: include 2 Pz orbitals. the symmetry is {I, RZ}x{I, M}. M is the mirror symmetry about the xy plane.
*/
enum class PointGroup {NONE,D6, C6, D4, D4m, D4m5, C4, D3, C3};
enum class ORBITAL {SINGLE, Dx2y2, Dz2, Px, Py, Pzu, Pzd};
enum class LINK_TYPE {SUPER_EXCHANGE_J, CHIRAL_K, HOPPING_T, CHARGE_TRANSFER_V, HUBBARD_U, EXCHANGE_J, PAIR_HOPPING_J};
enum class LADDER {PLUS, MINUS};
enum class SPIN {UP, DOWN};

/*
    ***************
    * Type Define *
    ***************
*/
// Matrix index and repI type
typedef unsigned long long int idx_t; 
typedef std::pair<idx_t,idx_t> pairIdx_t;
typedef unsigned int uint;
typedef std::complex<double> cdouble;
typedef cdouble dataType;


template<typename T>
using MAP = std::unordered_map<idx_t, T>;

typedef std::vector<int> VecI;
typedef std::vector<idx_t> VecIdx;
typedef std::vector<double> VecD;
typedef std::vector<cdouble> VecZ;

// link bond idx
typedef std::vector<std::vector<int>> BondMap;

/*
    ***************
    * Basis Class *
    ***************
*/
// use the binary form of an integer to represent basis state
#define BINARY_REP
#define KEEP_BASIS_NORM
// basis distributed among MPI workers
// #define DISTRIBUTED_BASIS
#ifdef DISTRIBUTED_BASIS
    // wavefunction vectors distributed among MPI workers
    #define DISTRIBUTED_STATE
#endif

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