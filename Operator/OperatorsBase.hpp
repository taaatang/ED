//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef OperatorsBase_hpp
#define OperatorsBase_hpp

#include "Operator/links.hpp"
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
    OperatorBase(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_ = 1, int dmNum_ = 0);
    virtual ~OperatorBase( ) { }
    
    OperatorBase& pushLink(Link<T> link, int matidx);
    OperatorBase& pushLinks(std::vector<Link<T>> links);
    void printLinks(bool brief = true) const;

    // Add onsite energy V
    virtual void pushV(std::vector<ORBITAL> orbList, double val) { }

    // Add onsite Coulomb interaction U
    virtual void pushU(std::vector<ORBITAL> orbList, double val) { }

    // O(t) --> O(t+dt)
    virtual bool next( ) {return false;}

    virtual void row(idx_t rowID, std::vector<MAP<T>>& rowMaps) = 0;
    
protected:
    Geometry *latt{nullptr};
    LATTICE_MODEL model{LATTICE_MODEL::HUBBARD};
    int linkCount{0};
    int spmCount{0};
    std::vector<Link<T>> Links, NCLinks;
    std::vector<Link<T>> superExchangeJ, chiralTermK;
    std::vector<Link<T>> hoppingT, interBandU, exchangeJ, pairHoppingJ;
};

template<typename T>
OperatorBase<T>::OperatorBase (Geometry *latt, Basis *Bi, Basis *Bf, int spmNum_, int dmNum_) :\
 FermionOperator<T>(Bi), SpinOperator<T>(Bi), SparseMatrix<T>(Bi, Bf, spmNum_, dmNum_){
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
#endif // OperatorsBase_hpp
