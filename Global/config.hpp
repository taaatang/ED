#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cmath>
#include <memory>
#include <algorithm>
#include "Global/globalPara.hpp"
#include "Utils/paras.hpp"
#include "Utils/path.hpp"
#include "Geometry/Geometry.hpp"

// std::string configFile = "../Input/config.txt";
// Parameters para(configFile);
std::string inputDir = "heisInput";

Parameters pathPara(inputDir, {"path.txt"});
Parameters para(inputDir, {"lattice.txt", "hamiltonian.txt"});
Parameters pulsePara(inputDir, {"pulse.txt"});
Parameters measurePara(inputDir, {"measure.txt"});

Path path(&pathPara, &para, &pulsePara);

std::unique_ptr<Geometry> latt;

std::unique_ptr<Basis> Bi, Bf;

std::unique_ptr<HamiltonianBase<dataType>> H, Hf;

template <typename T>
int sort(std::vector<T>& evals, std::vector<T*>& evecs, bool groundState = true, double degeneracyTol = 1e-8) {
    assert_msg(evals.size() == evecs.size(), "not have same number of eval and evec");
    using eigpair = std::pair<T,T*>;
    std::vector<eigpair> pairs;
    for (size_t i = 0; i < evals.size(); ++i) {
        pairs.push_back(std::make_pair(evals[i], evecs[i]));
    }
    std::sort(pairs.begin(), pairs.end(), [](const eigpair &a, const eigpair &b) {return std::real(a.first) < std::real(b.first);});
    for (size_t i = 0; i < pairs.size(); ++i) {
        evals[i] = pairs[i].first;
        evecs[i] = pairs[i].second;
    }
    int stateNum = 0;
    if (groundState) {
        auto w0 = std::real(evals.at(0));
        for (auto &w : evals) {
            if (std::abs(std::real(w) - w0) > degeneracyTol) break;
            ++stateNum;
        }
    } else {
        stateNum = evals.size();
    }
    return stateNum;
}


#endif // __CONFIG_H__