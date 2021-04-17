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
void sort(std::vector<T>& w, std::vector<T*>& state, bool groundState = true, double degeneracyTol = 1e-8) {
    using eigpair = std::pair<T,T*>;
    std::vector<eigpair> ws;
    for (size_t i = 0; i < w.size(); ++i) {
        ws.push_back(std::make_pair(w[i], state[i]));
    }
    std::sort(ws.begin(), ws.end(), [](const eigpair &a, const eigpair &b) {return std::real(a.first) < std::real(b.first);});
    if (groundState) {
        w.clear();
        state.clear();
        auto w0 = std::real(ws.at(0).first);
        for (auto &ep : ws) {
            if (std::abs(ep.first.real() - w0) > degeneracyTol) break;
            w.push_back(ep.first);
            state.push_back(ep.second);
        }
    } else {
        for (size_t i = 0; i < ws.size(); ++i) {
            w[i] = ws[i].first;
            state[i] = ws[i].second;
        }
    }
}


#endif // __CONFIG_H__