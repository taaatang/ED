#pragma once

#include <vector>
#include <string>
#include <optional>
#include <algorithm>
#include <concepts>

#include "basisState.hpp"
#include "geometry/geometry.hpp"
#include "utils/runtimeCheck.hpp"

template<typename T>
concept ElPhConfigEnabled = requires (T, int n, int nu, int nd, int nPho) { T::configure(n, nu, nd, nPho); };

template<typename T>
concept IsBasisType = std::derived_from<T, BasisStateInterface> && ElPhConfigEnabled<T>;

template<typename T>
concept IsPureSpin = std::same_as<T, SpinBasis>;

template<typename T>
concept IsPureCharge = std::same_as<T, ElectronBasis>;

template<typename T>
concept ContainPhonon =  std::same_as<T, ElectronPhononBasis>;

template<typename T>
concept ContainCharge = IsPureCharge<T> || std::same_as<T, ElectronPhononBasis>;

template<typename T>
concept ContainSpin = IsPureSpin<T> || ContainCharge<T>;

template<IsBasisType BasisState_t>
class Basis {
public:
    Basis(Geometry* lattptr, int nu, int nd, int maxPhPerSite, int k = -1, int p = -1);

    BasisState_t get(idx_t idx) const { return basisList.at(idx); }

    double norm(idx_t idx) const;

    std::optional<idx_t> search(const BasisState_t& state);

    void construct();

    void construct(int workerID, int workerNum);

    void construct(const std::string& basisFile, const std::string& normFile);

    void construct(const std::string& basisFile, const std::string& normFile, int workerID, int workerNum);

    bool empty() const;

    void clear();

    void print(std::ostream& os = std::cout, bool brief = true);

    template<IsBasisType T>
    friend std::ostream& operator<<(std::ostream& os, const Basis<T>& b);

private:
    bool isMinRep(BasisState_t state, double& basisNorm, const Generator<cdouble>* gs);

private:
    // total hilbert space dimension
	idx_t totDim{0};

    // lattice symmetry subspace dimension
	idx_t subDim{0};

	// local stored basis state number
	idx_t locDim{0};

    Geometry* latt{nullptr};

	int kIdx{-1};

	int pIdx{-1};

	std::vector<BasisState_t> basisList;

	std::vector<double> normList;

public:
    idx_t getTotDim() const { return totDim; }

    idx_t getSubDim() const { return subDim; }

    idx_t getLocDim() const { return locDim; }

    Geometry *getLatt() const { return latt; }

    int getKIdx() const { return kIdx; }

    void setKIdx(int kIdx_) { Basis::kIdx = kIdx_; }

    int getPIdx() const { return pIdx; }

    void setPIdx(int pIdx_) { Basis::pIdx = pIdx_; }

    bool isUseSymm() const { return kIdx != -1 || pIdx != -1; }
};

template<IsBasisType BasisState_t>
Basis<BasisState_t>::Basis(Geometry* lattptr, int nu, int nd, int maxPhPerSite, int k, int p) : latt(lattptr), kIdx(k), pIdx(p) {
    assert_msg(latt, "Null lattice_ptr!");
    BasisState_t::configure(latt->getOrbNum(), nu, nd, maxPhPerSite);
}

template<IsBasisType BasisState_t>
double Basis<BasisState_t>::norm(idx_t idx) const {
    if (isUseSymm()) {
        return normList.at(idx);
    } else {
        return 1.0;
    }
}

template<IsBasisType BasisState_t>
std::optional<idx_t> Basis<BasisState_t>::search(const BasisState_t &state) {
    auto low = std::lower_bound(basisList.begin(), basisList.end(), state);
    if (low != basisList.end() && *low == state) {
        return std::make_optional(low - basisList.begin());
    }
    return std::nullopt;
}

template<IsBasisType BasisState_t>
void Basis<BasisState_t>::construct() {
    clear();
    auto gs = latt->getTranslationGenerator(kIdx) * latt->getPointGroupGenerator(pIdx);
    BasisState_t state;
    totDim = state.getTotDim();
    if (isUseSymm()) {
        double basisNorm;
        do {
            if (isMinRep(state, basisNorm, &gs)) {
                basisList.push_back(state);
                normList.push_back(basisNorm);
            }
        } while (state.next());
    } else {
        do {
            basisList.push_back(state);
        } while (state.next());
        assert_msg(basisList.size() == totDim, "Wrong basis dimension!");
    }
    subDim = basisList.size();
    locDim = subDim;
}

//TODO
template<IsBasisType BasisState_t>
void Basis<BasisState_t>::construct(const std::string &basisFile, const std::string &normFile) {

}

//TODO
template<IsBasisType BasisState_t>
void
Basis<BasisState_t>::construct(const std::string &basisFile, const std::string &normFile, int workerID, int workerNum) {

}

//TODO:for hubbard model, we can use an empty basisList to save memory
template<IsBasisType BasisState_t>
bool Basis<BasisState_t>::empty() const {
    return basisList.empty();
}

template<IsBasisType BasisState_t>
bool Basis<BasisState_t>::isMinRep(BasisState_t state, double &basisNorm, const Generator<cdouble>* gs) {
    if (isUseSymm()) {
        cdouble cNorm2 = 0.0;
        for (const auto &u : gs->U) {
            BasisState_t trState = state;
            auto sgn = trState.transform(u);
            if (trState < state) {
                return false;
            } else if (trState == state) {
                cNorm2 += sgn * u.factor;
            }
        }
        // cNorm2 *= double(gs->G);
        assert_msg(std::imag(cNorm2) < INFINITESIMAL, "basis norm should ne real!");
        basisNorm = std::sqrt(std::real(cNorm2));
        if (basisNorm < INFINITESIMAL) {
            return false;
        }
        return true;
    } else {
        basisNorm = 1.0;
        return true;
    }
}

template<IsBasisType BasisState_t>
void Basis<BasisState_t>::clear() {
    basisList.clear();
    normList.clear();
}

template<IsBasisType BasisState_t>
void Basis<BasisState_t>::print(std::ostream& os, bool brief) {
    os << "k = " << getKIdx() << ", p = " << getPIdx() << ", total/sub dimension: "  << getTotDim() << "/" << getSubDim() << " = " << double(getTotDim())/getSubDim() << std::endl;
    if (!brief) {
        idx_t count = 0;
        for (const auto& state : basisList) {
            os << state << " norm = " << norm(count++) << std::endl;
        }
    }
}
