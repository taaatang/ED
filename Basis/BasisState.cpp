//
// Created by tatang on 8/1/21.
//

#include <bitset>
#include "basisState.hpp"
#include "utils/bitop.hpp"
#include "utils/io.hpp"
#include "utils/comb.hpp"

bool BasisStateInterface::configured = false;

int BasisStateInterface::nSite = 0;

int BasisStateInterface::nu = 0;

int BasisStateInterface::nd = 0;

int BasisStateInterface::maxPhPerSite = 0;

void BasisStateInterface::configureBase(int nSite_, int nu_, int nd_, int maxPhPerSite_) {
    if (!configured) {
        assert_msg(nSite_ >= 0 && nu_ >= 0 && nd_ >= 0 && maxPhPerSite_ >= 0, "nSite, nu, nd, maxPhPerSite should be non-negative!");
        assert_msg(nu_ <= nSite_ && nd_ <= nSite_, "nu and nd can't be greater than nSite!");
        nSite = nSite_;
        nu = nu_;
        nd = nd_;
        maxPhPerSite = maxPhPerSite_;
        configured = true;
    } else {
        std::cout << "Warning @ Attempt to reconfigure BasisStateInterface: Denied!" << std::endl;
    }
}

bool SpinBasis::isMinMaxSet = false;

BinaryState SpinBasis::min{};

BinaryState SpinBasis::max{};

SpinBasis::SpinBasis() {
    assert_msg(isConfigured() && isMinMaxSet, "SpinBasis is not configured!");
    state = min;
}

double SpinBasis::transform(const Transform<cdouble> &u) {
    spinState() = u.tr(spinState());
    return 1.0;
}

idx_t SpinBasis::getTotDim() {
    return BinaryState::getDim(getNSite(), getNu());
}

bool SpinBasis::next() {
    if (state < max) {
        state.next();
        return true;
    }
    return false;
}

void SpinBasis::configure(int nSite_, int nu_, int nd_, int maxPhPerSite_) {
    if (!isConfigured()) {
        BasisStateInterface::configureBase(nSite_, nu_, nSite_ - nu_, 0);
    }
    if (!isMinMaxSet) {
        min = BinaryState::min(getNSite(), getNu());
        max = BinaryState::max(getNSite(), getNu());
        isMinMaxSet = true;
    }
}

bool ElectronBasis::isMinMaxSet = false;

BinaryState ElectronBasis::upMin{};

BinaryState ElectronBasis::upMax{};

BinaryState ElectronBasis::dnMin{};

BinaryState ElectronBasis::dnMax{};

ElectronBasis ElectronBasis::min(BinaryState{}, BinaryState{});

ElectronBasis ElectronBasis::max(BinaryState{}, BinaryState{});

bool ElectronBasis::allowDoubleOcc = true;

ElectronBasis::ElectronBasis() {
    assert_msg(isConfigured() && isMinMaxSet, "ElectronBasis not configured!");
    *this = min;
    if (!allowDoubleOcc && isDoubleOcc()) {
        assert_msg(next(), "no valid state without double occupancy!");
    }
}

void ElectronBasis::configure(int nSite_, int nu_, int nd_, int maxPhPerSite_) {
    if (!allowDoubleOcc) {
        assert_msg(nu_ + nd_ <= nSite_, "nu + nd should not be greater than nSite when double occupancy is not allowed!");
    }
    if (!isConfigured()) {
        BasisStateInterface::configureBase(nSite_, nu_, nd_, 0);
    }
    if (!isMinMaxSet) {
        upMin = BinaryState::min(getNSite(), getNu());
        upMax = BinaryState::max(getNSite(), getNu());
        dnMin = BinaryState::min(getNSite(), getNd());
        dnMax = BinaryState::max(getNSite(), getNd());
        min = ElectronBasis(upMin, dnMin);
        max = ElectronBasis(upMax, dnMax);
        isMinMaxSet = true;
    }
}

bool ElectronBasis::next() {
    while (*this < max) {
        if (dn < dnMax) {
            dn.next();
        } else {
            dn = dnMin;
            up.next();
        }
//        std::cout << *this << std::endl;
        if (allowDoubleOcc) {
            return true;
        } else if (!isDoubleOcc()) {
            return true;
        }
    }
    return false;
}

double ElectronBasis::transform(const Transform<cdouble> &u) {
    double sgnUp;
    upState() = u.tr(upState(), sgnUp);
    double sgnDn;
    dnState() = u.tr(dnState(), sgnDn);
    return sgnUp * sgnDn;
}

idx_t ElectronBasis::getTotDim() {
    if (allowDoubleOcc) {
        return BinaryState::getDim(getNSite(), getNu()) * BinaryState::getDim(getNSite(), getNd());
    } else {
        assert_msg(nu + nd <= nSite, "invalid input for getTotDim");
        return BinaryState::getDim(nSite, nu) * BinaryState::getDim(nSite - nu, nd);
    }
}

bool ElectronPhononBasis::isMinMaxSet = false;

PhononState ElectronPhononBasis::phMin = PhononState{};

PhononState ElectronPhononBasis::phMax = PhononState{};

ElectronPhononBasis ElectronPhononBasis::min = ElectronPhononBasis(ElectronBasis{BinaryState{}, BinaryState{}}, PhononState{});

ElectronPhononBasis ElectronPhononBasis::max = ElectronPhononBasis(ElectronBasis{BinaryState{}, BinaryState{}}, PhononState{});

ElectronPhononBasis::ElectronPhononBasis() {
    assert_msg(isConfigured() && isMinMaxSet, "ElectronPhononBasis is not configured!");
    el = ElectronBasis();
    ph = phMin;
}

void ElectronPhononBasis::configure(int nSite_, int nu_, int nd_, int maxPhPerSite_) {
    if(!isConfigured()) {
        configureBase(nSite_, nu_, nd_, maxPhPerSite_);
    }
    ElectronBasis::configure(nSite_, nu_, nd_, maxPhPerSite_);
    if (!isMinMaxSet) {
        phMin = PhononState::min(getNSite());
        phMax = PhononState::max(getNSite(), getMaxPhPerSite());
        min = ElectronPhononBasis(ElectronBasis::getMin(), phMin);
        max = ElectronPhononBasis(ElectronBasis::getMax(), phMax);
        isMinMaxSet = true;
    }
}

bool ElectronPhononBasis::next() {
    if (*this < max) {
        if (ph < phMax) {
            ph.next(maxPhPerSite);
            return true;
        } else {
            if (el.next()) {
                ph = phMin;
                return true;
            }
            return false;
        }
    }
    return false;
}

double ElectronPhononBasis::transform(const Transform<cdouble> &u) {
    ph() = u.tr(ph());
    return el.transform(u);
}

idx_t ElectronPhononBasis::getTotDim() {
    return PhononState::getDim(getNSite(), maxPhPerSite) * el.getTotDim();
}

BinaryState::BinaryState(int nSite, int nParticle) {
    *this = BinaryState::min(nSite, nParticle);
}

void BinaryState::next() {
    state = nextLexicographicalNumber(state);
}

BinaryState BinaryState::min(int nSite, int nParticle) {
    BinaryState b;
    for (int i = 0; i < nParticle; ++i) {
        bitSet(b(), i);
    }
    return b;
}

BinaryState BinaryState::max(int nSite, int nParticle) {
    BinaryState b;
    for (int i = 0; i < nParticle; ++i) {
        bitSet(b(), nSite - 1 - i);
    }
    return b;
}

idx_t BinaryState::getDim(int nSite, int nParticle) {
    return combination(idx_t(nSite), idx_t(nParticle));
}

void BinaryState::print(std::ostream &os, int n) const {
    for (--n; n >= 0; --n) {
        if (count(n)) {
            os << 1;
        } else {
            os << 0;
        }
    }
}

PhononState::PhononState(int nSite, uint8_t maxPhononPerSite) {
    *this = PhononState::min(nSite);
}

void PhononState::next(int maxPhPerSite) {
    for (auto rit = state.rbegin(); rit != state.rend(); ++rit) {
        if (*rit < maxPhPerSite) {
            *rit += 1;
            break;
        } else {
            *rit = 0;
        }
    }
}

PhononState PhononState::min(int nSite) {
    PhononState ph;
    ph.state = std::vector<uint8_t>(nSite, 0);
    return ph;
}

PhononState PhononState::max(int nSite, int nMaxPhPerSite) {
    PhononState ph;
    ph.state = std::vector<uint8_t>(nSite, nMaxPhPerSite);
    return ph;
}

idx_t PhononState::getDim(int nSite, int nMaxPhPerSite) {
    idx_t dim = 1;
    idx_t dimPerSite = nMaxPhPerSite + 1;
    for (int i = 0; i < nSite; ++i) {
        dim *= dimPerSite;
    }
    return dim;
}

bool operator<(const BinaryState &lhs, const BinaryState &rhs) {
    return lhs.state < rhs.state;
}

bool operator==(const BinaryState &lhs, const BinaryState &rhs) {
    return lhs.state == rhs.state;
}

bool operator<(const PhononState &lhs, const PhononState &rhs) {
    assert_msg(lhs.state.size() == rhs.state.size(), "mismatch basis state length!");
    for (size_t i = 0; i < lhs.state.size(); ++i) {
        if (lhs.state[i] < rhs.state[i]) {
            return true;
        } else if (lhs.state[i] > rhs.state[i]) {
            return false;
        }
    }
    return false;
}

bool operator==(const PhononState &lhs, const PhononState &rhs) {
    assert_msg(lhs.state.size() == rhs.state.size(), "mismatch basis state length!");
    for (size_t i = 0; i < lhs.state.size(); ++i) {
        if (lhs.state[i] != rhs.state[i]) {
            return false;
        }
    }
    return true;
}

bool operator<(const SpinBasis &lhs, const SpinBasis &rhs) {
    return lhs.state < rhs.state;
}

bool operator==(const SpinBasis &lhs, const SpinBasis &rhs) {
    return lhs.state == rhs.state;
}

bool operator<(const ElectronBasis &lhs, const ElectronBasis &rhs) {
    if ((lhs.up < rhs.up) || (lhs.up == rhs.up && lhs.dn < rhs.dn)) {
        return true;
    }
    return false;
}

bool operator==(const ElectronBasis &lhs, const ElectronBasis &rhs) {
    return lhs.up == rhs.up && lhs.dn == rhs.dn;
}

bool operator<(const ElectronPhononBasis &lhs, const ElectronPhononBasis &rhs) {
    if ((lhs.el < rhs.el) || (lhs.el == rhs.el && lhs.ph < rhs.ph)) {
        return true;
    }
    return false;
}

bool operator==(const ElectronPhononBasis &lhs, const ElectronPhononBasis &rhs) {
    return lhs.el == rhs.el && lhs.ph == rhs.ph;
}

std::ostream &operator<<(std::ostream &os, BinaryState b) {
    os << std::bitset<64>(b.state);
    return os;
}

std::ostream& operator<<(std::ostream &os, const PhononState &p) {
    os << "ph:";
    for (auto n : p.state) {
        os << int(n);
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const SpinBasis &s) {
    os << "spin:";
    s.state.print(os, s.getNSite());
    return os;
}

std::ostream &operator<<(std::ostream &os, const ElectronBasis &e) {
    os << "up:";
    e.up.print(os, e.getNSite());
    os << " dn:";
    e.dn.print(os, e.getNSite());
    return os;
}

std::ostream &operator<<(std::ostream &os, const ElectronPhononBasis &ep) {
    os << ep.el;
    os << ' ' << ep.ph;
    return os;
}
