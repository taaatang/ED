#pragma once

#include <vector>
#include <complex>
#include <utility> // pair
#include <unordered_map>
#include <string>

using idx_t = uint64_t;
using pairIdx_t = std::pair<idx_t,idx_t>;
using cdouble = std::complex<double>;
using dataType = cdouble;

using VecI = std::vector<int>;
using VecIdx = std::vector<idx_t>;
using VecD = std::vector<double>;
using VecZ = std::vector<cdouble>;
using ArrD = std::vector<VecD>;

using Str = std::string;
using VecStr = std::vector<std::string>;

/*
    ***************
    * Basis Class *
    ***************
*/
// basis distributed among MPI workers
// #define DISTRIBUTED_BASIS