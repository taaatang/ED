#ifndef __GLOBALTYPE_H__
#define __GLOBALTYPE_H__

#include <vector>
#include <complex>
#include <utility> // pair
#include <unordered_map>
#include <string>

enum class LATTICE_MODEL {HUBBARD, tJ, HEISENBERG, ElPh};
/*
    D4m: for multi-band hubbard, due to orbital phases, the symmetry is only a sub group {I, RZ}. R is rotation, Z is reflection.
    D4m5: include 2 Pz orbitals. the symmetry is {I, RZ}x{I, M}. M is the mirror symmetry about the xy plane.
*/
enum class PointGroup {NONE,D6, C6, D4, D4m, D4m5, C4, D3, C3};
enum class ORBITAL {SINGLE, Dx2y2, Dz2, Px, Py, Pzu, Pzd};
/**
 * @brief Interaction Type
 * 
 */
enum class LINK_TYPE {SZ, N, NUND, SUPER_EXCHANGE_J, CHIRAL_K, HOPPING_T, CHARGE_TRANSFER_V, HUBBARD_U, EXCHANGE_J, PAIR_HOPPING_J, NCHARGE_SITE_PHONON, NCHARGE_BOND_PHONON, HOPPING_BOND_PHONON};
enum class LADDER {PLUS, MINUS};
enum class SPIN {UP, DOWN};

/*
    ***************
    * Type Define *
    ***************
*/
// Matrix index and repI type
using idx_t =  unsigned long long int; 
using pairIdx_t = std::pair<idx_t,idx_t>;
using uint = unsigned int;
using cdouble = std::complex<double>;
using dataType = cdouble;


template<typename T>
using MAP = std::unordered_map<idx_t, T>;

using VecI = std::vector<int>;
using VecIdx = std::vector<idx_t>;
using VecD = std::vector<double>;
using VecZ = std::vector<cdouble>;
using ArrD = std::vector<VecD>;

using Str = std::string;
using VecStr = std::vector<std::string>;


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