#pragma once
//
// Created by tatang on 8/3/21.
//

#include <optional>
#include "global/typeAlias.hpp"
#include "basis/basis.hpp"
#include "utils/bitop.hpp"

enum class SPIN {UP, DOWN};

enum class LADDER {PLUS, MINUS};

template<typename T, IsBasisType B>
struct BasisVal {
    B basis;
    T val{1.0};
    BasisVal() = delete;
    BasisVal(B b) : basis(b) {}
};

template<typename T, IsBasisType B>
using BVopt = std::optional<BasisVal<T, B>>;

template<SPIN S>
struct CPlus {
    int pos{0};
    explicit CPlus(int i = 0) : pos(i) {}
    const CPlus& operator()(int i) {
        pos = i;
        return *this;
    }
};

template<SPIN S>
struct CMinus {
    int pos{0};
    explicit CMinus(int i = 0) : pos(i) { }
    const CMinus& operator()(int i) {
        pos = i;
        return *this;
    }
};

template<SPIN S>
struct NCharge {
    int pos{0};
    explicit NCharge(int i = 0) : pos(i) { }
    const NCharge& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct APlus {
    int pos{0};
    explicit APlus(int i = 0) : pos(i) { }
    const APlus& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct AMinus {
    int pos{0};
    explicit AMinus(int i = 0) : pos(i) { }
    const AMinus& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct NPhonon {
    int pos{0};
    explicit NPhonon(int i = 0) : pos(i) { }
    const NPhonon& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct SPlus {
    int pos{0};
    explicit SPlus(int i = 0) : pos(i) { }
    const SPlus& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct SMinus {
    int pos{0};
    explicit SMinus(int i = 0) : pos(i) { }
    const SMinus& operator()(int i) {
        pos = i;
        return *this;
    }
};

struct Sz {
    int pos{0};
    explicit Sz(int i = 0) : pos(i) { }
    const Sz& operator()(int i) {
        pos = i;
        return *this;
    }
};

template<SPIN S, ContainCharge B, typename T>
BVopt<T, B> operator*(const CPlus<S>& cp, BVopt<T, B> bv) {
    if (bv) {
        idx_t* state;
        if constexpr(S == SPIN::UP) {
            state = &bv->basis.upState();
        } else {
            state = &bv->basis.dnState();
        }
        if(!bitTest(*state, cp.pos)) {
            bitFlip(*state, cp.pos);
            int counter = bitCount(*state, cp.pos);
            if constexpr(S == SPIN::DOWN) {
                counter += bitCount(bv->basis.upState(), B::getNSite());
            }
            int sgn = (counter % 2 == 0) ? 1 : -1;
            bv->val *= sgn;
        } else {
            return std::nullopt;
        }
    }
    return bv;
}

template<SPIN S, ContainCharge B, typename T>
BVopt<T, B> operator*(const CMinus<S>& cm, BVopt<T, B> bv) {
    if (bv) {
        idx_t* state;
        if constexpr(S == SPIN::UP) {
            state = &bv->basis.upState();
        } else {
            state = &bv->basis.dnState();
        }
        if(bitTest(*state, cm.pos)) {
            bitFlip(*state, cm.pos);
            int counter = bitCount(*state, cm.pos);
            if constexpr(S == SPIN::DOWN) {
                counter += bitCount(bv->basis.upState(), B::getNSite());
            }
            int sgn = (counter % 2 == 0) ? 1 : -1;
            bv->val *= sgn;
        } else {
            return std::nullopt;
        }
    }
    return bv;
}

template<SPIN S, ContainCharge B, typename T>
BVopt<T, B> operator*(const NCharge<S>& nch, BVopt<T, B> bv) {
    if (bv) {
        idx_t* state;
        if constexpr(S == SPIN::UP) {
            state = &bv->basis.upState();
        } else {
            state = &bv->basis.dnState();
        }
        if(bitTest(*state, nch.pos)) {
            bv->val *= 1.0;
        } else {
            bv->val *= 0.0;
        }
    }
    return bv;
}

template<ContainPhonon B, typename T>
BVopt<T, B> operator*(const APlus& ap, BVopt<T, B> bv) {
    if (bv) {
        auto& state = bv->basis.phState();
        if (state.at(ap.pos) < B::getMaxPhPerSite()) {
            state[ap.pos] += 1;
            bv->val *= std::sqrt(double(state[ap.pos]));
            return bv;
        }
    }
    return std::nullopt;
}

template<ContainPhonon B, typename T>
BVopt<T, B> operator*(const AMinus& a, BVopt<T, B> bv) {
    if (bv) {
        auto& state = bv->basis.phState();
        if (state.at(a.pos) > 0) {
            bv->val *= std::sqrt(double(state[a.pos]));
            state[a.pos] -= 1;
            return bv;
        }
    }
    return std::nullopt;
}

template<ContainPhonon B, typename T>
BVopt<T, B> operator*(const NPhonon& n, BVopt<T, B> bv) {
    if (bv) {
        bv->val *= double(bv->basis.phState().at(n.pos));
        return bv;
    }
    return std::nullopt;
}

template<ContainSpin B, typename T>
BVopt<T, B> operator*(const SPlus& sp, BVopt<T, B> bv) {
    if (bv) {
        if constexpr (IsPureSpin<B>) {
            auto& state = bv->basis.spinState();
            if (!bitTest(state, sp.pos)) {
                bitFlip(state, sp.pos);
                return bv;
            }
        } else if constexpr (ContainCharge<B>){
            return CPlus<SPIN::UP>(sp.pos) * (CMinus<SPIN::DOWN>(sp.pos) * bv);
        }
    }
    return std::nullopt;
}

template<ContainSpin B, typename T>
BVopt<T, B> operator*(const SMinus& sm, BVopt<T, B> bv) {
    if (bv) {
        if constexpr (IsPureSpin<B>) {
            auto& state = bv->basis.spinState();
            if (bitTest(state, sm.pos)) {
                bitFlip(state, sm.pos);
                return bv;
            }
        } else if constexpr (ContainCharge<B>){
            return CPlus<SPIN::DOWN>(sm.pos) * (CMinus<SPIN::UP>(sm.pos) * bv);
        }
    }
    return std::nullopt;
}

template<ContainSpin B, typename T>
BVopt<T, B> operator*(const Sz& sz, BVopt<T, B> bv) {
    if (bv) {
        if constexpr (IsPureSpin<B>) {
            auto& state = bv->basis.spinState();
            bv->val *= bitTest(state, sz.pos) ? 0.5 : -0.5;
            return bv;
        } else if constexpr (ContainCharge<B>){
            auto& up = bv->basis.upState();
            auto& dn = bv->basis.dnState();
            if (bitTest(up, sz.pos) != bitTest(dn, sz.pos)) {
                bv->val *= bitTest(up, sz.pos) ? 0.5 : -0.5;
                return bv;
            }
        }
    }
    return std::nullopt;
}

template<typename T, IsBasisType B>
BVopt<T, B> operator*(T fac, BVopt<T, B> bv) {
    if (bv) {
        bv->val *= fac;
        return bv;
    }
    return std::nullopt;
}

inline std::ostream& operator<<(std::ostream& os, SPIN s) {
    // os<<"spin_";
    switch (s) {
        case SPIN::UP:
            os<<"up";
            break;
        case SPIN::DOWN:
            os<<"dn";
            break;
        default:
            os<<"undefined";
            break;
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& os, LADDER t) {
    // os<<"ladder_";
    switch (t) {
        case LADDER::PLUS:
            os<<"plus";
            break;
        case LADDER::MINUS:
            os<<"minus";
            break;
        default:
            os<<"undefined";
            break;
    }
    return os;
}