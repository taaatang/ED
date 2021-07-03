//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef OperatorsBase_hpp
#define OperatorsBase_hpp

#include <algorithm>
#include "Operator/links.hpp"
#include "Operator/interaction_transform.hpp"
#include "Operator/Operations.hpp"
#include "Operator/SparseMatrix.hpp"
#include "Global/globalPara.hpp"
#include "Basis/Basis.hpp"
#include "Geometry/Geometry.hpp"
#include "Pulse/pulse.hpp"

template<typename T>
class OperatorBase: public FermionOperator<T>, public SpinOperator<T>, public SparseMatrix<T>{
public:
    OperatorBase( ) { }
    OperatorBase(Geometry* latt, Basis* Bi, Basis* Bf, bool commuteWithTrans = false, bool commuteWithPG = false, int spmNum_ = 1, int dmNum_ = 0);
    virtual ~OperatorBase( ) { }
    
    OperatorBase& pushLink(Link<T> link, int matidx);
    OperatorBase& pushLinks(std::vector<Link<T>> links);
    void printLinks(bool brief = true) const;
    // generate symmetry transformed operators
    void getGiGf(Generator<T>& Gi, Generator<T>& Gf, std::vector<Transform<T>>& allTr) const;
    void transform();

    // Add onsite energy V
    virtual void pushV(std::vector<ORBITAL> orbList, double val) { }

    // Add onsite Coulomb interaction U
    virtual void pushU(std::vector<ORBITAL> orbList, double val) { }

    // O(t) --> O(t+dt)
    virtual bool next( ) {return false;}

    virtual void row(idx_t rowID, std::vector<MAP<T>>& rowMaps) = 0;
    
protected:
    Geometry *latt{nullptr};
    bool commuteWithTrans{false};
    bool commuteWithPG{false};
    LATTICE_MODEL model{LATTICE_MODEL::HUBBARD};
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

template<typename T>
OperatorBase<T>::OperatorBase (Geometry *latt, Basis *Bi, Basis *Bf, bool commuteWithTrans_, bool commuteWithPG_, int spmNum_, int dmNum_) :\
 FermionOperator<T>(Bi, commuteWithTrans_, commuteWithPG_), SpinOperator<T>(Bi, commuteWithTrans, commuteWithPG), SparseMatrix<T>(Bi, Bf, spmNum_, dmNum_){
    commuteWithTrans = commuteWithTrans_;
    commuteWithPG = commuteWithPG_;
    this->latt = latt;
    this->model = this->Bi->getModel();
    assert(this->Bi->getModel() == this->Bf->getModel());
}

template<typename T>
OperatorBase<T>& OperatorBase<T>::pushLink(Link<T> link, int matID){
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

template<typename T>
OperatorBase<T>& OperatorBase<T>::pushLinks(std::vector<Link<T>> links){
    assert(spmCount<SparseMatrix<T>::spmNum);
    for (auto& link : links) pushLink(link, spmCount);
    if (links.at(0).isOrdered()) spmCount += 2;
    else spmCount++;
    assert(spmCount<=SparseMatrix<T>::spmNum);
    return *this;
}

template<typename T>
void OperatorBase<T>::printLinks(bool brief) const {
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


template<typename T>
void OperatorBase<T>::getGiGf(Generator<T>& Gi, Generator<T>& Gf, std::vector<Transform<T>>& allTr) const {
    Gi.setIdentity(latt->getOrbNum());
    Gf = latt->getGP(this->Bf->getPGIndex()) * latt->getGT(this->Bf->getkIndex());
    if (commuteWithPG) {
        Gf = latt->getGP(this->Bi->getPGIndex()) * Gf;
    } else {
        Gi = latt->getGP(this->Bi->getPGIndex()) * Gi;
    }
    if (commuteWithTrans) {
        Gf = latt->getGT(this->Bi->getkIndex()) * Gf;
    } else {
        Gi = latt->getGT(this->Bi->getkIndex()) * Gi;
    }
    auto Gtot = Gi * Gf; 
    allTr = Gtot.U;
}

template<typename T>
void OperatorBase<T>::transform() {
    Generator<T> Gi, Gf;
    std::vector<Transform<T>> allTr;
    getGiGf(Gi, Gf, allTr);
    assignTrInteractions<T, 2>(Gi, Gf, allTr, superExchangeJ, trSuperExchangeJ, 'n');
    assignTrInteractions<T, 3>(Gi, Gf, allTr, chiralTermK, trChiralTermK, 'c');
}
#endif // OperatorsBase_hpp
