#pragma once

#include <vector>
#include <memory>
#include "operator/operators.hpp"
#include "utils/paras.hpp"
#include "utils/path.hpp"
#include "solver/spectra.hpp"
#include "solver/timeEvolver.hpp"
#include "utils/timer.hpp"

template <typename T, IsBasisType B>
void timeEvolve(Hamiltonian<T, B>* H, const T* sate, double dt, int steps, std::vector<OperatorBase<T, B>*> ops);

template <typename T, IsBasisType B>
void doubleTimeCorrelator(OperatorBase<T, B>* opt, OperatorBase<T, B>* opi, Hamiltonian<T, B>* H, std::vector<T>& sate, double dt, int steps, const std::string& dir);

template <typename T, IsBasisType B>
void dynamicLanczos();

template <typename T, IsBasisType B>
void dynamicBicg();


template <typename T, IsBasisType B>
void doubleTimeCorrelator(OperatorBase<T, B>* op, Hamiltonian<T, B>* Hi, Hamiltonian<T, B>* Hf, const T* state, int krydim, double dt, int steps, const std::string& dir, bool isMaster) {
    std::vector<T> results;
    std::vector<T> stateL(state, state + Hi->getnloc());
    std::vector<T> stateR(Hf->getnloc(), 0.0);
    op->MxV(stateL.data(), stateR.data());
    results.push_back(std::conj(op->vMv(stateR.data(), stateL.data())));
    TimeEvolver right(stateR.data(), Hf, krydim);
    TimeEvolver left(stateL.data(), Hi, krydim);
    int progressPrint = steps / 20;
    if (isMaster) std::cout << "Begin time corr: |";
    for (int count = 0; count < steps; ++count) {
        right.evolve(-dt);
        left.evolve(-dt);
        results.push_back(std::conj(op->vMv(right.getVec(), left.getVec())));
        if (isMaster) {
            if (count % progressPrint == 0) {
                std::cout << ">" << std::flush;
            }
        }
    }
    if (isMaster) std::cout << "done!" << std::endl;
    if (isMaster) {
        mkdir_fs(dir);
        save<double>(&dt, 1, dir + "/dt");
        save<T>(results.data(), int(results.size()), dir + "/corr");
    }
}
