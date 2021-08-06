//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//
#pragma once

#include <algorithm>
#include "Global/globalPara.hpp"
#include "Geometry/Geometry.hpp"
#include "Operator/SparseMatrix.hpp"
#include "Operator/links.hpp"
#include "Operator/interaction_transform.hpp"
#include "Operator/localOperators.hpp"
#include "Pulse/pulse.hpp"

// push data to an unordered map
template <typename T>
inline void MapPush(MAP<T>* map_pt, idx_t key, T val){
    auto it = map_pt->find(key);
    if constexpr (std::is_same<cdouble, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = std::conj(val);
        } else {
            it->second += std::conj(val);
        }
    } else if constexpr (std::is_same<double, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = val;
        } else {
            it->second += val;
        }
    } else {
        std::cout<<"MapPush only defined for double/cdouble!\n";
        exit(1);
    }
}

template<typename T, IsBasisType B>
class OperatorBase: public SparseMatrix<T> {
public:
    OperatorBase( ) = delete;

    OperatorBase(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool commuteWithTrans = false, bool commuteWithPG = false, int spmNum_ = 1, int dmNum_ = 0);

    virtual ~OperatorBase( ) = default;
    
    OperatorBase& pushLink(Link<T> link, int matidx);

    OperatorBase& pushLinks(std::vector<Link<T>> links);

    void printLinks(bool brief = true) const;

    // generate symmetry transformed operators
    void getGiGf(Generator<T>& Gi, Generator<T>& Gf, std::vector<Transform<T>>& allTr) const;

    OperatorBase& transform();

    // Add onsite energy V
    virtual OperatorBase& pushV(std::vector<ORBITAL> orbList, double val) { return *this; }

    // Add onsite Coulomb interaction U
    virtual OperatorBase& pushU(std::vector<ORBITAL> orbList, double val) { return *this; }

    // O(t) --> O(t+dt)
    virtual bool next( ) {return false;}

    void pushElement(const BVopt<T, B>& bv, MAP<T>* rowMap);

    virtual void row(idx_t rowID, std::vector<MAP<T>>& rowMaps) = 0;
    
protected:
    Geometry *latt{nullptr};
    Basis<B>* Bi{nullptr};
    Basis<B>* Bf{nullptr};
    bool commuteWithTrans{false};
    bool commuteWithPG{false};
    int linkCount{0};
    int spmCount{0};
    std::vector<Link<T>> Links, NCLinks;
    // Heisenberg terms
    std::vector<Link<T>> superExchangeJ, chiralTermK;
    Interactions<T,1> Sz{LINK_TYPE::SZ};
    std::vector<TrInteractions<T,1>> trSz;
    std::vector<TrInteractions<T,2>> trSuperExchangeJ;
    std::vector<TrInteractions<T,3>> trChiralTermK;
    // Hubbard terms
    std::vector<Link<T>> hoppingT, interBandU, exchangeJ, pairHoppingJ;
    std::vector<TrInteractions<T,2>> trHoppingT, trInterBandU, trExchangeJ, trPairHoppingJ;
};

template<typename T, IsBasisType B>
OperatorBase<T, B>::OperatorBase (Geometry *latt_, Basis<B> *Bi_, Basis<B> *Bf_, bool commuteWithTrans_, bool commuteWithPG_, int spmNum_, int dmNum_) : SparseMatrix<T>(Bf_->getSubDim(), Bi_->getSubDim(), spmNum_, dmNum_),
        latt(latt_), Bi(Bi_), Bf(Bf_), commuteWithTrans(commuteWithTrans_), commuteWithPG(commuteWithPG_) {

}

template<typename T, IsBasisType B>
OperatorBase<T, B>& OperatorBase<T, B>::pushLink(Link<T> link, int matID){
    if(matID==0)assert(link.isConst());
    else assert(!link.isConst());
    link.genLinkMaps(latt);
    if (link.getLinkNum() > 0) {
        link.setid(linkCount, matID);
        switch (link.getLinkType()) {
            case LINK_TYPE::SUPER_EXCHANGE_J:
                superExchangeJ.push_back(link);
                break;
            case LINK_TYPE::CHIRAL_K:
                chiralTermK.push_back(link);
                break;
            case LINK_TYPE::HOPPING_T:
                hoppingT.push_back(link);
                break;
            case LINK_TYPE::HUBBARD_U:
                interBandU.push_back(link);
                break;
            case LINK_TYPE::EXCHANGE_J:
                exchangeJ.push_back(link);
                break;
            case LINK_TYPE::PAIR_HOPPING_J:
                pairHoppingJ.push_back(link);
                break;
            default:
                std::cout<<"link_type not defined!\n";
                exit(1);
                break;
        }
        linkCount++;
    }
    return *this;
    // if(matID==0)assert(link.isConst());
    // else assert(!link.isConst());
    // Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(latt); 
    // linkCount++;
    // return *this;
}

template<typename T, IsBasisType B>
OperatorBase<T, B>& OperatorBase<T, B>::pushLinks(std::vector<Link<T>> links){
    assert(spmCount<SparseMatrix<T>::spmNum);
    for (auto& link : links) pushLink(link, spmCount);
    if (links.at(0).isOrdered()) spmCount += 2;
    else spmCount++;
    assert(spmCount<=SparseMatrix<T>::spmNum);
    return *this;
}

template<typename T, IsBasisType B>
void OperatorBase<T, B>::printLinks(bool brief) const {
    for (const auto& link:hoppingT) {
        link.print(brief);
    }
    for (const auto& link:superExchangeJ) {
        link.print(brief);
    }
    for (const auto& link:chiralTermK) {
        link.print(brief);
    }
}

template<typename T, size_t N>
void addTrInteractions(const Generator<T>& Gi, const Generator<T>& Gf, const Link<T>& link, std::vector<TrInteractions<T, N>>& trOp) {
    for (size_t i = 0; i < Gi.G; ++i) {
        Generator<T> Gs;
        Gs.add(Gi[i]);
        auto op = linkToTrInteractions<T, N>(Gi[i], link); 
        auto Gtot = Gs * Gf;
        for (size_t j = 0; j < Gtot.G; ++j) {
            TrInteractions<T, N> comp(Gtot[j], link.getLinkType());
            auto found = std::lower_bound(trOp.begin(), trOp.end(), comp);
            if (*found == comp) {
                found->Op += op * Gtot[j].factor;
            } else {
                assert_msg(false, "transformation not found in addTrOp!");
            }
        }
    }
}

template<typename T, size_t N>
void assignTrInteractions(const Generator<T>& Gi, const Generator<T>& Gf, const std::vector<Transform<T>>& totTr, const std::vector<Link<T>>& links, std::vector<TrInteractions<T, N>>& trOps, char order) {
    trOps.clear();
    if (links.empty()) {
        return;
    }
    for (const auto& t : totTr) {
        trOps.push_back(TrInteractions<T,N>(t, links[0].getLinkType()));
    }
    for (const auto& link : links) {
        addTrInteractions<T, N>(Gi, Gf, link, trOps);
    }
    for(auto& trOp : trOps) {
        trOp.Op.condense(order);
    }
}

template<typename T, size_t N>
void addTrInteractions(const Generator<T>& Gi, const Generator<T>& Gf, const Interactions<T, N>& Op, std::vector<TrInteractions<T, N>>& trOps) {
    for (size_t i = 0; i < Gi.G; ++i) {
        Generator<T> Gs;
        Gs.add(Gi[i]);
        auto trOp = Gi[i] * Op; 
        auto Gtot = Gs * Gf;
        for (size_t j = 0; j < Gtot.G; ++j) {
            TrInteractions<T, N> comp(Gtot[j], Op.type);
            auto found = std::lower_bound(trOps.begin(), trOps.end(), comp);
            if (*found == comp) {
                found->Op += trOp * Gtot[j].factor;
            } else {
                assert_msg(false, "transformation not found in addTrOp!");
            }
        }
    }
}

template<typename T, size_t N>
void assignTrInteractions(const Generator<T>& Gi, const Generator<T>& Gf, const std::vector<Transform<T>>& totTr, const std::vector<Interactions<T,N>>& Ops, std::vector<TrInteractions<T, N>>& trOps, char order) {
    trOps.clear();
    if (Ops.empty()) {
        return;
    }
    for (const auto& t : totTr) {
        trOps.push_back(TrInteractions<T,N>(t, Ops[0].type));
    }
    for (const auto& Op : Ops) {
        addTrInteractions<T, N>(Gi, Gf, Op, trOps);
    }
    for(auto& trOp : trOps) {
        trOp.Op.condense(order);
    }
}


template<typename T, IsBasisType B>
void OperatorBase<T, B>::getGiGf(Generator<T>& Gi, Generator<T>& Gf, std::vector<Transform<T>>& allTr) const {
    Gi.setIdentity(latt->getOrbNum());
    Gf = latt->getPointGroupGenerator(this->Bf->getPIdx()) * latt->getTranslationGenerator(this->Bf->getKIdx());
    if (commuteWithPG) {
        Gf = latt->getPointGroupGenerator(this->Bi->getPIdx()) * Gf;
    } else {
        Gi = latt->getPointGroupGenerator(this->Bi->getPIdx()) * Gi;
    }
    if (commuteWithTrans) {
        Gf = latt->getTranslationGenerator(this->Bi->getKIdx()) * Gf;
    } else {
        Gi = latt->getTranslationGenerator(this->Bi->getKIdx()) * Gi;
    }
    auto Gtot = Gi * Gf; 
    allTr = Gtot.U;
}

template<typename T, IsBasisType B>
OperatorBase<T, B>& OperatorBase<T, B>::transform() {
    Generator<T> Gi, Gf;
    std::vector<Transform<T>> allTr;
    getGiGf(Gi, Gf, allTr);
    assignTrInteractions<T, 2>(Gi, Gf, allTr, superExchangeJ, trSuperExchangeJ, 'n');
    assignTrInteractions<T, 3>(Gi, Gf, allTr, chiralTermK, trChiralTermK, 'c');
    for (const auto& trOp : trSuperExchangeJ) {
        std::cout << trOp.g.factor << " * " << trOp.g.transformList << std::endl;
        for (const auto& bond : trOp.Op.bonds) {
            std::cout << bond.val << " * " << bond.sites[0] << "--" << bond.sites[1] << std::endl;
        }
    }
    return *this;
}

template<typename T, IsBasisType B>
void OperatorBase<T, B>::pushElement(const BVopt<T, B> &bv, MAP<T> *rowMap) {
    if (bv) {
        auto idxOpt = Bi->search(bv->basis);
        if (idxOpt) {
            auto ni = Bi->norm(*idxOpt);
            MapPush(rowMap, *idxOpt, bv->val / ni);
        }
    }
}
