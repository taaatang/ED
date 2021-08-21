#pragma once

#include <vector>
#include <complex>
#include <utility> // pair
#include <unordered_map>
#include <string>

/**
 * @brief Interaction Type
 * 
 */
enum class LINK_TYPE {SZ, N, NUND, SUPER_EXCHANGE_J, CHIRAL_K, HOPPING_T, CHARGE_TRANSFER_V, HUBBARD_U, EXCHANGE_J, PAIR_HOPPING_J, PHONON_W0, NCHARGE_SITE_PHONON, NCHARGE_BOND_PHONON, HOPPING_BOND_PHONON};


/*
    ***************
    * Type Define *
    ***************
*/
// Matrix index and repI type
using idx_t = uint64_t;
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
// basis distributed among MPI workers
// #define DISTRIBUTED_BASIS

/*
    *******************
    * LANCZOSIterator *
    *******************
*/

enum LANCZOS_OPTION {ALPHA_BETA, ALPHA_BETA_Q};
