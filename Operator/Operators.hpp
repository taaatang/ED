//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Operators_hpp
#define Operators_hpp

#include "../Operator/OperatorsBase.hpp"
#include "../Operator/SparseMatrix.hpp"

/* Operators_hpp */


/*
    ***********
    * Hubbard *
    ***********
*/

template <class T>
class Hubbard: public FermionOperator, public SparseMatrix<T>{
private:
    // constant links
    int linkCount;
    int spmCount;
    std::vector<Link<T>> Links;
    // non-constant links
    std::vector<Link<T>> NCLinks;
    // V:charge transfer energy on different orbitals. U:onsite repulsion
    VecD V, U;
    Geometry *pt_lattice;
public:
    Hubbard(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int dmNum_=1, int spindim=2):V(pt_lat->getUnitOrbNum()), U(pt_lat->getUnitOrbNum()),\
        FermionOperator(pt_Ba),SparseMatrix<T>(pt_Ba, pt_Ba, pt_Ba->getSubDim(),spmNum_,dmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
    ~Hubbard(){};

    Hubbard& pushLink(Link<T> link, int matID){
        if(matID==0)assert(link.isConst());
        else assert(!link.isConst());
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); 
        linkCount++;
        return *this;
    }
    Hubbard& pushLinks(std::vector<Link<T>> Links_){
        assert(spmCount<SparseMatrix<T>::spmNum);
        for (int i = 0; i < Links_.size(); i++) pushLink(Links_[i], spmCount);
        if (Links_.at(0).isOrdered()) spmCount += 2;
        else spmCount++;
        assert(spmCount<=SparseMatrix<T>::spmNum);
        return *this;
    }
    Hubbard& pushV(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++){VecI ids = pt_lattice->getOrbID(*it); for(auto id=ids.begin();id!=ids.end();id++)V.at(*id)=val;} return *this;}
    Hubbard& pushU(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++){VecI ids = pt_lattice->getOrbID(*it); for(auto id=ids.begin();id!=ids.end();id++)U.at(*id)=val;} return *this;}
    void printV() const {std::cout<<"V:"<<std::endl;for(auto it=V.begin();it!=V.end();it++)std::cout<<*it<<", ";std::cout<<std::endl;}
    void printU() const {std::cout<<"U:"<<std::endl;for(auto it=U.begin();it!=U.end();it++)std::cout<<*it<<", ";std::cout<<std::endl;}
    double diagVal(VecI occ, VecI docc) const {double val=0.0; for(int i=0;i<pt_lattice->getUnitOrbNum();i++)val += occ.at(i)*V.at(i)+docc.at(i)*U.at(i);return val;}
    void setVal(int matID, T val){(SparseMatrix<T>::parameters).at(matID) = val;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    void genMat();
};

class Current: public FermionOperator, public SparseMatrix<cdouble>{
    int linkCount;
    int spmCount;
    std::vector<Link<cdouble>> Links;
    Geometry *pt_lattice;
    std::string plz;
public:
    Current(Geometry *pt_lat, Basis *pt_Ba, std::string plzIn, int spmNum_=1):FermionOperator(pt_Ba),SparseMatrix<cdouble>(pt_Ba, pt_Ba, pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0),plz(plzIn){}
    ~Current(){};
    Current& pushLink(Link<cdouble> link, int matID){
        // int idx = 0; if(plz=="x") idx = 0; else if(plz=="y") idx = 1; else if(plz=="z") idx=2; else exit(1);
        // auto linkVec = link.getvec(0);
        // double sign = 1.0; if(linkVec.at(idx)>0.) sign = -1.0; else if(linkVec.at(idx)==0.) sign = 0.0;
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); //Links[linkCount].setVal(link.getVal()*sign);
        linkCount++;
        return *this;
    }
    Current& pushLinks(std::vector<Link<cdouble>> Links_){
        for (int i = 0; i < Links_.size(); i++) pushLink(Links_[i], 0);
        spmCount++;
        assert(spmCount<=SparseMatrix<cdouble>::spmNum);
        return *this;
    }
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
};

class Nocc: public SparseMatrix<cdouble>{
private:
    Geometry* pt_lattice;
    Basis* pt_Basis;
public:
    Nocc(Geometry *pt_lat, Basis *pt_Ba);
    ~Nocc(){}

    void row(ind_int rowID, std::vector<MAP>& rowMaps){};
    inline void row(ind_int rowID);
    void genMat();
    double count(ORBITAL orbital, dataType* vec);
};

/*
    ******
    * tJ *
    ******
*/

template <class T>
class HtJ: public FermionOperator, public SpinOperator, public SparseMatrix<T>{
private:
    // constant links
    int linkCount;
    int spmCount;
    std::vector<Link<T>> Links;
    // non-constant links
    std::vector<Link<T>> NCLinks;
    Geometry *pt_lattice;
    Basis* pt_Basis;
public:
    HtJ(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int dmNum_=0, int spindim=2):\
        FermionOperator(pt_Ba),\
        SpinOperator(pt_Ba),\
        SparseMatrix<T>(pt_Ba, pt_Ba, pt_Ba->getSubDim(),spmNum_,dmNum_), pt_lattice(pt_lat), pt_Basis(pt_Ba),linkCount(0),spmCount(0){}
    ~HtJ(){}

    HtJ& pushLink(Link<T> link, int matID){
        if(matID==0)assert(link.isConst());
        else assert(!link.isConst());
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); 
        linkCount++;
        return *this;
    }
    HtJ& pushLinks(std::vector<Link<T>> Links_){
        assert(spmCount<SparseMatrix<T>::spmNum);
        for (int i = 0; i < Links_.size(); i++) pushLink(Links_[i], spmCount);
        if (Links_.at(0).isOrdered()) spmCount += 2;
        else spmCount++;
        assert(spmCount<=SparseMatrix<T>::spmNum);
        return *this;
    }
    void setVal(int matID, T val){(SparseMatrix<T>::parameters).at(matID) = val;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    void genMat();
};

/*
    **************
    * Heisenberg *
    **************
*/

template <class T>
class Heisenberg: public SpinOperator, public SparseMatrix<T>{
private:
    // constant links
    int linkCount;
    int spmCount;
    std::vector<Link<T>> Links;
    // non-constant links
    std::vector<Link<T>> NCLinks;
    Geometry *pt_lattice;
public:
    Heisenberg(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2):\
        SpinOperator(pt_Ba),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
    ~Heisenberg(){};

    Heisenberg& pushLink(Link<T> link, int matID){
        if(matID==0)assert(link.isConst());
        else assert(!link.isConst());
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); 
        linkCount++;
        return *this;
    }
    Heisenberg& pushLinks(std::vector<Link<T>> Links_){
        assert(spmCount<SparseMatrix<T>::spmNum);
        for (int i = 0; i < Links_.size(); i++) pushLink(Links_[i], spmCount);
        spmCount++;
        assert(spmCount<=SparseMatrix<T>::spmNum);
        return *this;
    }

    void setVal(int matID, double val){(SparseMatrix<T>::parameters).at(matID) = val;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    void genMat();
};

template<class T>
class RamanOp: public SpinOperator, public SparseMatrix<T>{
private:
    int linkCount;
    int spmCount;
    std::vector<Link<T>> Links;
    std::vector<Link<T>> NCLinks;
    Geometry *pt_lattice;

    // polarization of in/out photons
    VecD pIn, pOut;
    std::vector<VecD> RamanWeight;
    bool NoRW;
public:
    RamanOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2):\
        SpinOperator(pt_Ba),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){
            pIn = VecD{1.0,0.0,0.0};
            pOut = VecD{0.0,1.0,0.0};
            NoRW = true;        
        }
    ~RamanOp(){};

    RamanOp& pushLink(Link<T> link, int matID){
        if(matID==0)assert(link.isConst());
        else assert(!link.isConst());
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); 
        linkCount++;
        return *this;
    }
    RamanOp& pushLinks(std::vector<Link<T>> Links_){
        assert(spmCount<SparseMatrix<T>::spmNum);
        for (int i = 0; i < Links_.size(); i++) pushLink(Links_[i], spmCount);
        spmCount++;
        assert(spmCount<=SparseMatrix<T>::spmNum);
        return *this;
    }

    void setplz(VecD pIn_, VecD pOut_);
    void setVal(int matID, double val){(SparseMatrix<T>::parameters).at(matID) = val;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
};


template <class T>
class SzkOp: public SpinOperator, public SparseMatrix<cdouble>{
private:
    // initial k index Ki and final k index Kf
    int Ki, Kf;
    Basis *pt_Bi,*pt_Bf;
    Geometry *pt_lattice;
    std::vector<cdouble> expFactor;
public:
    SzkOp(Geometry *pt_lat, Basis *pt_Bi_, Basis *pt_Bf_, int spmNum_=1, int spindim=2):pt_Bi(pt_Bi_),pt_Bf(pt_Bf_), pt_lattice(pt_lat),expFactor(pt_lattice->getSiteNum()),\
        SpinOperator(pt_Bi_),SparseMatrix<cdouble>(pt_Bi_,pt_Bf_,pt_Bf_->getSubDim(),spmNum_){
            assert(pt_Bi->getPGIndex()==-1 and pt_Bf->getPGIndex()==-1);
            Ki = pt_Bi->getkIndex();
            Kf = pt_Bf->getkIndex();
            // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
            for (int i = 0; i < pt_lattice->getSiteNum(); i++) {
                expFactor[i] = pt_lattice->expKR(Kf,i)/pt_lattice->expKR(Ki,i);
            }
        };
    ~SzkOp(){};
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    // void genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q);
    void genMat();
};

// correlator:Si*Si+r
// can also be used to construct total spin operator: S^2
template <class T>
class SSOp: public SpinOperator, public SparseMatrix<T>{
private:
    Geometry *pt_lattice;
    int r; // if r = -1, sum.r sum.i Si*Si+r. total S^2
    std::vector<VecI> siteJList;
public:
    SSOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2);
    ~SSOp(){}
    void setr(int r_){assert(r_<pt_lattice->getSiteNum());r=r_;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    // S(i)*S(i+r)
    void genPairMat(int rIndex);
    void project(double s, T* vec);
};

/*
    ****************************
    * Operator Implementations *
    ****************************
*/

template <class T>
void Hubbard<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    // diagonal part. occupancy and double-occ
    VecI occ, docc;
    ind_int repI = pt_Basis->getRepI(rowID);
    pairIndex pairRepI = pt_Basis->getPairRepI(repI);
    pt_lattice->orbOCC(pairRepI, occ, docc);
    T val = diagVal(occ,docc);
    for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
        if((*linkit).getLinkType()!=LINK_TYPE::HUBBARD_U)continue;
        for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
            int siteI = (*bondit).at(0);
            int siteJ = (*bondit).at(1);
            assert(siteI!=siteJ);
            val += (*linkit).getVal()*(bitTest(pairRepI.first,siteI)+bitTest(pairRepI.second,siteI))*(bitTest(pairRepI.first,siteJ)+bitTest(pairRepI.second,siteJ));
        }
    }
    SparseMatrix<T>::putDiag(val,rowID);
    // off diagonal part
    std::vector<ind_int> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); i++){
        pairIndex pairRepI = pt_Basis->getPairRepI(finalIndList[i]);
        bool isfminRep = pt_Basis->isfMin(pairRepI.first);
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            if((*linkit).getLinkType()!=LINK_TYPE::HOPPING_T)continue;
            int matID = (*linkit).getmatid();
            int matIDp = matID; 
            if ((*linkit).isOrdered()) matIDp++;
            cdouble factor = factorList.at(i) * (*linkit).getVal();
            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                int siteI = (*bondit).at(0);
                int siteJ = (*bondit).at(1);
                // cp.siteI * cm.siteJ
                cpcm(SPIN::SPIN_UP, siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                cpcm(SPIN::SPIN_UP, siteJ, siteI, factor, finalIndList[i], &rowMaps[matIDp]);
                if(isfminRep){
                    cpcm(SPIN::SPIN_DOWN, siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                    cpcm(SPIN::SPIN_DOWN, siteJ, siteI, factor, finalIndList[i], &rowMaps[matIDp]);   
                }
            }
        }
    }
    // diag(rowID,val,&rowMaps[0]);
}

template <class T>
void Hubbard<T>::genMat(){
    int kIndex = pt_Basis->getkIndex();
    SparseMatrix<T>::clear();
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < spmNum; i++) SparseMatrix<T>::pushRow(&rowMap,i);
    VecI initVec(pt_lattice->getOrbNum()), initVecp(pt_lattice->getOrbNum());
    for (ind_int rowID = BaseMatrix<T>::startRow; rowID < BaseMatrix<T>::endRow; rowID++){
        /*
            *****************
            * Constant Part *
            * ***************
        */
        rowMap.clear();
        // diagonal part. occupancy and double-occ
        VecI occ, docc;
        ind_int initInd = pt_Basis->getRepI(rowID);
        pt_Basis->repToVec(initInd, initVec, initVecp);
        pt_lattice->orbOCC(initVec, initVecp, occ, docc);
        double val = diagVal(occ,docc);
        cdouble diag_val = 0.0;
        // off diagonal part
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->repToVec(finalIndList[i], initVec, initVecp);
            if(finalIndList[i]==initInd)diag_val += val*factorList[i];
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                    int siteI = (*bondit).at(0);
                    int siteJ = (*bondit).at(1);
                    // cp.siteI * cm.siteJ
                    cpcm(SPIN::SPIN_UP, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_UP, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMap);   
                }
            }
        }
        diag(rowID,diag_val,&rowMap);
        SparseMatrix<T>::pushRow(&rowMap);

        /*
            ***********************
            * Time Dependent Part *
            * *********************
        */
        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->repToVec(finalIndList[i], initVec, initVecp);
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // cm.siteID * cp.siteIDp
                    cpcm(SPIN::SPIN_UP, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_UP, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap); 
                }
            }
            SparseMatrix<T>::pushRow(&rowMap, (*linkit).getmatid());
        }
    }
}
inline void Nocc::row(ind_int rowID){
    // diagonal part. occupancy
    VecI occ;
    ind_int repI = pt_Basis->getRepI(rowID);
    pairIndex pairRepI = pt_Basis->getPairRepI(repI);
    pt_lattice->orbOCC(pairRepI, occ);
    ind_int loc_rowID = rowID - BaseMatrix<cdouble>::startRow;
    for (int i = 0; i < pt_lattice->getUnitOrbNum(); i++) {SparseMatrix<cdouble>::diagValList[i][loc_rowID] = occ[i];}
}

template <class T>
void HtJ<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    // off diagonal part
    std::vector<ind_int> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); i++){
        pairIndex pairRepI = pt_Basis->getPairRepI(finalIndList[i]);
        bool isfminRep = pt_Basis->isfMin(pairRepI.first);
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            LINK_TYPE type = (*linkit).getLinkType();
            int matID = (*linkit).getmatid();
            int matIDp = matID; 
            if ((*linkit).isOrdered()) matIDp++;
            cdouble factor = factorList.at(i) * (*linkit).getVal();
            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                int siteI = (*bondit).at(0);
                int siteJ = (*bondit).at(1);
                switch(type){
                    case LINK_TYPE::HOPPING_T:{
                        // cp.siteI * cm.siteJ
                        cpcm(SPIN::SPIN_UP, siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                        cpcm(SPIN::SPIN_UP, siteJ, siteI, factor, finalIndList[i], &rowMaps[matIDp]);
                        if(isfminRep){
                            cpcm(SPIN::SPIN_DOWN, siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                            cpcm(SPIN::SPIN_DOWN, siteJ, siteI, factor, finalIndList[i], &rowMaps[matIDp]);   
                        }
                        break;
                    }
                    case LINK_TYPE::SUPER_EXCHANGE_J:{
                        // szi * szj - 1/4*ni*nj 
                        szsznn(siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                        // 1/2 * sm.siteID * sp.siteIDp
                        spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[matID]);
                        // 1/2 * sp.siteID * sm.siteIDp
                        smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[matID]);
                        break;
                    }
                    default:{
                        std::cout<<"Interaction type not defined for tJ (HtJ::row)\n";
                        exit(1);
                    }
                }   
            }
        }
    }
    // diag(rowID,val,&rowMaps[0]);
}

template <class T>
void RamanOp<T>::setplz(VecD pIn_, VecD pOut_){
    pIn = pIn_;
    pOut = pOut_;
    RamanWeight.resize(Links.size());
    double norm = std::sqrt((pt_lattice->RdotR(pIn,pIn))*(pt_lattice->RdotR(pOut,pOut)));
    for(int linkid=0; linkid<Links.size(); linkid++){
        RamanWeight[linkid].resize(Links[linkid].getLinkVecNum());
        for(int vecid=0; vecid<Links[linkid].getLinkVecNum(); vecid++){
            VecD r = Links[linkid].getvec(vecid);
            RamanWeight[linkid][vecid] = (pt_lattice->RdotR(pIn,r))*(pt_lattice->RdotR(pOut,r))/(pt_lattice->RdotR(r,r))/norm;
        }
    }
    NoRW = false;
}
template <class T>
void RamanOp<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    if(NoRW) setplz(pIn, pOut);
    std::vector<ind_int> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); i++){
        int linkid = 0;
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            int matID = (*linkit).getmatid();
            int matIDp = matID; 
            if ((*linkit).isOrdered()) matIDp++;
            cdouble factor = factorList.at(i) * (*linkit).getVal();
            int bondid = 0;
            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                double rw = RamanWeight.at(linkid).at((*linkit).getvecid(bondid));
                cdouble factor_rw = factor * rw;
                int siteI = (*bondit).at(0);
                int siteJ = (*bondit).at(1);
                // szi * szj - 1/4*ni*nj 
                // szsznn(siteI, siteJ, factor_rw, finalIndList[i], &rowMaps[matID]);
                szsz(siteI, siteJ, factor_rw, finalIndList[i], &rowMaps[matID]);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteI, siteJ, factor_rw/2.0, finalIndList[i], &rowMaps[matID]);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteI, siteJ, factor_rw/2.0, finalIndList[i], &rowMaps[matID]);
                bondid++;
            }
            linkid++;
        }
    }
}

template <class T>
void Heisenberg<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    if(pt_Basis->getSiteDim()==2){
        int kIndex = pt_Basis->getkIndex();
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                int matID = (*linkit).getmatid();
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                switch((*linkit).getLinkType()){
                    case LINK_TYPE::SUPER_EXCHANGE_J:{
                        for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                            int siteID = (*bondit).at(0);
                            int siteIDp = (*bondit).at(1);
                            // sz.siteID * sz.siteIDp
                            szsz(siteID, siteIDp, factor, finalIndList[i], &rowMaps[matID]);
                            // 1/2 * sm.siteID * sp.siteIDp
                            spsm(siteID, siteIDp, factor/2.0, finalIndList[i], &rowMaps[matID]);
                            // 1/2 * sp.siteID * sm.siteIDp
                            smsp(siteID, siteIDp, factor/2.0, finalIndList[i], &rowMaps[matID]);
                        }
                        break;
                    }

                    case LINK_TYPE::CHIRAL_K:{
                        for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                            int siteI = (*bondit).at(0);
                            int siteJ = (*bondit).at(1);
                            int siteK = (*bondit).at(2);
                            chiral(siteI, siteJ, siteK, factor, finalIndList[i], &rowMaps[matID]);
                        }
                        break;
                    }

                    default: break;
                }
                
            }
        }
    }
    else{
        int kIndex = pt_Basis->getkIndex();
        VecI initVec(pt_lattice->getOrbNum());
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->repToVec(finalIndList[i], initVec);
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                int matID = (*linkit).getmatid();
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMaps[matID]);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
                }
            }
        }
    }  
}
// generate Hamiltonian in the subspacd labeled by kIndex
template <class T>
void Heisenberg<T>::genMat(){
    int kIndex = pt_Basis->getkIndex();
    SparseMatrix<T>::clear();
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < SparseMatrix<T>::spmNum; i++) SparseMatrix<T>::pushRow(&rowMap,i);
    VecI initVec(pt_lattice->getOrbNum());
    // calculate <R1k|H*Pk|R2k>/norm1/norm2
    for (ind_int rowID = SparseMatrix<T>::startRow; rowID < SparseMatrix<T>::endRow; rowID++){
        rowMap.clear();
        std::vector<cdouble> factorList;
        std::vector<ind_int> finalIndList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->repToVec(finalIndList[i], initVec);
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                }
            }
        }
        SparseMatrix<T>::pushRow(&rowMap);

        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->repToVec(finalIndList[i], initVec);
                cdouble factor = factorList.at(i) * (*linkit).getVal();
                for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                }
            }
            SparseMatrix<T>::pushRow(&rowMap, (*linkit)->getmatid());
        }
    }
}

/*
    ********************
    * Correlator Class *
    ********************
*/
template <class T>
SSOp<T>::SSOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_, int spindim):r(-1),pt_lattice(pt_lat), siteJList(pt_lat->getSiteNum()),\
    SpinOperator(pt_Ba),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_){
    VecD coordi(3), coordr(3), coordf(3);
    for (int rIndex = 0; rIndex < pt_lat->getSiteNum();rIndex++){
        siteJList.at(rIndex).resize(pt_lat->getSiteNum());
        pt_lattice->getSiteR(rIndex, coordr.data());
        for (int i = 0; i < pt_lattice->getOrbNum(); i++){
            pt_lattice->getOrbR(i,coordi.data());
            vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), 3);
            int siteJ;
            if (pt_lattice->coordToOrbid(coordf.data(), siteJ)){
                siteJList.at(rIndex).at(i) = siteJ;
            }else{
                std::cout<<"translation position not found for orbid = "<<i<<", transVecid = "<<rIndex<<std::endl;
                exit(1);
            }
        }
    }
}

template <class T>
void SSOp<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    //binary rep
    if(pt_Basis->getSiteDim()==2){
        int kIndex = pt_Basis->getkIndex();
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); siteI++){
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    szsz(siteI, siteJ, factor, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                }
                else{
                    for (int rIndex = 0; rIndex < pt_lattice->getSiteNum(); rIndex++){
                        int siteJ = siteJList[rIndex][siteI];
                        // sz.siteI * sz.siteJ
                        szsz(siteI, siteJ, factor, finalIndList[i], &rowMaps[0]);
                        // 1/2 * sm.siteI * sp.siteJ
                        spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                        // 1/2 * sp.siteI * sm.siteJ
                        smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                    }
                }
            }
        }
    }
    else{
        int kIndex = pt_Basis->getkIndex();
        VecI initVec(pt_lattice->getOrbNum());
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->repToVec(finalIndList[i], initVec);
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); siteI++){
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                }
                else{
                    for (int rIndex = 0; rIndex < pt_lattice->getSiteNum(); rIndex++){
                        int siteJ = siteJList[rIndex][siteI];
                        // sz.siteI * sz.siteJ
                        szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMaps[0]);
                        // 1/2 * sm.siteI * sp.siteJ
                        spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                        // 1/2 * sp.siteI * sm.siteJ
                        smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                    }
                }
            }
        }
    }
}
template <class T>
void SSOp<T>::project(double s, T* vec){
    assert(r==-1);
    double stot = s*(s+1);
    std::vector<T> tmp(BaseMatrix<T>::get_nloc());
    double smin = double(pt_lattice->getOrbNum()%2)/2.0;
    double smax = double(pt_lattice->getOrbNum())/2.0+0.1;
    for(double si = smin; si<smax; si += 1.0){
        if (std::abs(si-s)>0.1){
            MxV(vec, tmp.data());
            double stoti = si * (si+1.0);
            #pragma omp parallel for
            for(ind_int i =0; i < BaseMatrix<T>::get_nloc();i++) vec[i] = (tmp[i]-stoti*vec[i])/(stot-stoti);
        }
    }
}
// generate matrix in subsapce labeled by kIndex for sum.r:Sr*Sr+dr, dr is labeled by rIndex
template <class T>
void SSOp<T>::genPairMat(int rIndex){
    #ifdef DISTRIBUTED_BASIS
        exit(1);
    #else
    assert(rIndex>0 and rIndex<pt_lattice->getSiteNum());
    r = rIndex;
    int kIndex = pt_Basis->getkIndex();
    clear();
    reserve(pt_lattice->getOrbNum()/2+1);
    MAP rowMap;
    pushRow(&rowMap);
    VecI initVec(pt_lattice->getOrbNum());
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->repToVec(finalIndList[i], initVec);
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); siteI++){
                int siteJ = siteJList[rIndex][siteI];
                // sz.siteI * sz.siteJ
                szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMap);
                // 1/2 * sm.siteI * sp.siteJ
                spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
                // 1/2 * sp.siteI * sm.siteJ
                smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
            }
        }
        pushRow(&rowMap);
    }
    #endif
}

template <class T>
void SzkOp<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    if(pt_Bi->getSiteDim()==2){
        #ifdef DISTRIBUTED_BASIS
            ind_int repI = pt_Bf->getRepI(rowID);
            cdouble dval = 0.0;
            for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                dval += getSz(siteID,repI) * expFactor[siteID];
            }
            dval *= pt_Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            ind_int colID;
            if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                // std::cout<<"colID:"<<colID<<"\n";
                ind_int repI = pt_Bi->getRepI(colID);
                cdouble dval = 0.0;
                for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                    dval += getSz(siteID,repI) * expFactor[siteID];
                }
                dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
                // std::cout<<"rowMaps size:"<<rowMaps.size()<<"\n";
                rowMaps[0][colID] = dval;
            }
        #endif
    }
    else{
        #ifdef DISTRIBUTED_BASIS
            ind_int repI = pt_Bf->getRepI(rowID);
            VecI initVec(pt_lattice->getOrbNum());
            cdouble dval = 0.0;
            pt_Bf->repToVec(repI, initVec);
            for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                dval += getSz(siteID,initVec) * expFactor[siteID];
            }
            dval *= pt_Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            ind_int colID;
            if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                VecI initVec(pt_lattice->getOrbNum());
                cdouble dval = 0.0;
                pt_Bi->repToVec(pt_Bi->getRepI(colID), initVec);
                for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                    dval += getSz(siteID,initVec) * expFactor[siteID];
                }
                dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
                rowMaps[0][colID] = dval;
            }
        #endif
    }
}
template <class T>
void SzkOp<T>::genMat(){
    #ifdef DISTRIBUTED_BASIS
        exit(1);
    #else
    clear();
    reserve(1);
    MAP rowMap;
    pushRow(&rowMap);
    cdouble dval;
    VecI initVec(pt_lattice->getOrbNum());
    switch(PARTITION){
        case ROW_PARTITION:{
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                rowMap.clear();
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bi->repToVec(pt_Bi->getRepI(colID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += getSz(siteID,initVec) * expFactor[siteID];
                    }
                    dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
                    rowMap[colID] = dval;
                }
                pushRow(&rowMap);
            }
            break;
        }
        // col Partition need to be checked!
        case COL_PARTITION:{
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bf->repToVec(pt_Bf->getRepI(rowID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += getSz(siteID,initVec) * expFactor[siteID];
                    }
                    dval *= pt_Bi->getNorm(rowID)/pt_Bf->getNorm(colID);
                    rowMap[colID] = dval;
                }
                pushRow(&rowMap);
            }
            break;
        }
    }
    #endif
}

#endif // Operators_hpp
