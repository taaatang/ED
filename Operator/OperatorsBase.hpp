//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef OperatorsBase_hpp
#define OperatorsBase_hpp

#include "OperatorBasics.hpp"
#include "SparseMatrix.hpp"
#include "../Global/globalPara.hpp"
#include "../Basis/Basis.hpp"
#include "../Geometry/Geometry.hpp"

template<typename T>
class OperatorBase: public FermionOperator<T>, public SpinOperator<T>, public SparseMatrix<T>{
public:
    OperatorBase( ) { }
    OperatorBase(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_ = 1, int dmNum_ = 0);
    virtual ~OperatorBase( ) { }
    
    OperatorBase& pushLink(Link<T> link, int matidx);
    OperatorBase& pushLinks(std::vector<T> links);
    virtual void row(idx_t rowID, std::vector<MAP>& rowMaps) = 0;
protected:
    Basis *Bi{nullptr}, *Bf{nullptr};
    Geometry *latt{nullptr};
    LATTICE_MODEL model{HUBBARD};
    int linkCount{0};
    int spmCount{0};
    std::vector<Link<T>> Links, NCLinks;
};

template<typename T>
OperatorBase<T>::Operator(Geometry *latt, Basis *Bi, Basis *Bf, int spmNum_, int dmNum_) :\
 FermionOperator<T>(Bi), SpinOperator<T>(Bi), SparseMatrix<T>(Bi, Bf, Bf->getSubDim(), spmNum_. dmNum_){
    this->latt = latt;
    this->model = Bi->getModel();
    assert(model == MODEL && Bi->getModel() == Bf->getModel());
}

template<typename T>
OperatorBase<T>& OperatorBase<T>::pushLink(Link<T> link, int matID){
    if(matID==0)assert(link.isConst());
    else assert(!link.isConst());
    Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(latt); 
    linkCount++;
    return *this;
}

template<typename T>
OperatorBase<T>& OperatorBase<T>::pushLinks(std::vector<T> links){
    assert(spmCount<SparseMatrix<T>::spmNum);
    for (auto& link : links) pushLink(link, spmCount);
    if (links.at(0).isOrdered()) spmCount += 2;
    else spmCount++;
    assert(spmCount<=SparseMatrix<T>::spmNum);
    return *this;
}

#endif // OperatorsBase_hpp
