//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Operators_hpp
#define Operators_hpp

#include "../Basis/Basis.hpp"
#include "SparseMatrix.hpp"

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>

/* Operators_hpp */


template <class T>
class Link{
    LINK_TYPE link_type;
    T val;
    // same group ids will be contained in the same sparse matrix
    int linkid,matid;
    bool is_Const, is_Ordered;
    T (*timeFunc)(int);
    bool is_timeFunc_set;
    int timeStep;
    std::vector<ORBITAL> orbList;
    std::vector<VecD> LinkVecs;
    std::vector<VecI> LinkList;
public:
    Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, bool is_Const_ = true, bool is_Ordered_ = false):timeStep(0), is_timeFunc_set(false){
        link_type = link_type_;
        orbList = orbList_;
        val = val_;
        is_Const = is_Const_;
        is_Ordered = is_Ordered_;
    };
    ~Link(){};
    
    LINK_TYPE getLinkType() const {return link_type;}
    // link value
    T getVal() const {return val;}
    int getlinkid() const {return linkid;}
    int getmatid() const {return matid;}
    int getLinkNum() const {return LinkList.size();}
    // rep orbital of the link
    ORBITAL getRepOrb() const {return orbList.at(0);}
    // number of orbs in one link
    int getLinkSize() const {return orbList.size()-1;}
    // size of LinkVecs
    int getLinkVecNum() const {return LinkVecs.size();}
    // is this a const link or time dependent link
    bool isConst() const {return is_Const;}
    bool isOrdered() const {return is_Ordered;}
    // begin/end iterator for LinkList
    std::vector<VecI>::iterator begin() {return LinkList.begin();}
    std::vector<VecI>::iterator end() {return LinkList.end();}

    void setVal(T val_){val = val_;}
    void setid(int linkid_, int matid_=0){linkid = linkid_; matid = matid_;}
    // set time evolve function for val
    void setTimeFunc(T (*timeFunc_)(int)){timeFunc = timeFunc_; is_timeFunc_set = true;}
    void setVal(){
        if(!is_Const) {
            if(!is_timeFunc_set){std::cout<<"time evolve function not set for linkid:"<<linkid<<std::endl; exit(1);}
            val = timeFunc(timeStep); 
        }
        timeStep++;
    }
    // add one link (contains orbids in the same link)
    void push_back(std::vector<int> orbidList){LinkList.push_back(orbidList);}
    // add a orb's relative coord to the rep orb
    Link<T>& addLinkVec(std::vector<double> vec){assert(vec.size()==3); LinkVecs.push_back(vec); return *this;}
    //generate all the links->linkList
    void genLinkMaps(Geometry* pt_lattice);
    void print() const;
};

template <class T>
void Link<T>::genLinkMaps(Geometry* pt_lattice){
    VecD coordi(pt_lattice->getDim()),coordf(pt_lattice->getDim()); 
    for (int orbid = 0; orbid < pt_lattice->getOrbNum(); orbid++){
        if(pt_lattice->is_Orbital(orbid,getRepOrb())){
            pt_lattice->getOrbR(orbid,coordi.data());
            for (int j = 0; j < getLinkVecNum(); j+=getLinkSize()){
                VecI tmp;
                tmp.push_back(orbid);
                bool is_bond = true;
                for (int k = j; k < j+getLinkSize(); k++){
                    vecXAdd(1.0,coordi.data(),1.0,LinkVecs.at(k).data(),coordf.data(),3);
                    int orbidf;
                    if (pt_lattice->coordToOrbid(coordf.data(), orbidf)){
                        assert(pt_lattice->getOrb(orbidf) == orbList.at(k-j+1));
                        tmp.push_back(orbidf);
                    }else{
                        is_bond = false;
                        break;
                    }
                }
                if(is_bond)push_back(tmp);
            }
        }   
    }
}

template <class T>
void Link<T>::print() const {
    std::cout<<"Link id:"<<linkid<<", value:"<<val<<std::endl;
    int counter = 0;
    for (auto bondit = LinkList.begin(); bondit != LinkList.end(); bondit++){
        counter++;
        std::cout<<"Bond "<<counter<<": ";
        for(auto it = (*bondit).begin(); it != (*bondit).end(); it++){
            std::cout<<*it<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

/*
    *******************
    * Operators Class *
    *******************
*/
class FermionOperator{
protected:
    Basis* pt_Basis;
    LATTICE_MODEL model;
public:
    FermionOperator(Basis* pt_Ba, LATTICE_MODEL mod=LATTICE_MODEL::HUBBARD):pt_Basis(pt_Ba),model(mod){};
    ~FermionOperator(){};
    void diag(ind_int rowID, dataType factor, MAP* rowMap){
        #ifdef DISTRIBUTED_BASIS
            ind_int repI = pt_Basis->getRepI(rowID);
            MapPush(rowMap,repI,factor);
        #else
            MapPush(rowMap,rowID,factor);
        #endif
    }

    bool cpcm(int siteI, int siteJ, ind_int repI, ind_int repIp, ind_int& repIf, int& sign){
        repIf = repI;
        if (siteI==siteJ){
            if (bitTest(repI,siteJ)){
                sign = 1;
                return true;
            }
        }else{
            switch(model){
                case HUBBARD:{
                    if(bitTest(repI,siteJ) && (!bitTest(repI,siteI))){
                        bitFlip(repIf,siteI);
                        bitFlip(repIf,siteJ);
                        int counter = 0;
                        if (siteI < siteJ) {
                            for (int i = siteI + 1; i < siteJ; i++) if (bitTest(repI,i)) counter++;
                        }else{
                            for (int i = siteJ + 1; i < siteI; i++) if (bitTest(repI,i)) counter++;
                        }
                        if (counter%2==0) sign = 1;
                        else sign = -1;
                        return true;
                    }
                    break;
                }

                // no double occupancy
                case t_J:{
                    if(bitTest(repI,siteJ) && (!bitTest(repI,siteI)) && (!bitTest(repIp,siteI))){
                        bitFlip(repIf,siteI);
                        bitFlip(repIf,siteJ);
                        int counter = 0;
                        if (siteI < siteJ) {
                            for (int i = siteI + 1; i < siteJ; i++) if (bitTest(repI,i)) counter++;
                        }else{
                            for (int i = siteJ + 1; i < siteI; i++) if (bitTest(repI,i)) counter++;
                        }
                        if (counter%2==0) sign = 1;
                        else sign = -1;
                        return true;
                    }
                    break;
                }
                default:break;
            }
        }
        return false;
    }
    void cpcm(SPIN spin, int siteI, int siteJ, dataType factor, ind_int repI, MAP* rowMap){
        pairIndex pairRepI=pt_Basis->getPairRepI(repI);
        pairIndex pairRepIf = pairRepI;
        int sign;
        ind_int colidx;
        switch(spin){
            case SPIN_UP:{
                if (cpcm(siteI, siteJ, pairRepI.first, pairRepI.second, pairRepIf.first, sign)){
                    #ifdef DISTRIBUTED_BASIS
                    if(pt_Basis->isfMin(pairRepIf.first)) MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                    #else
                    if (pt_Basis->search(pairRepIf, colidx)){
                        dataType val = factor * sign;
                        double finalNorm = pt_Basis->getNorm(colidx);
                        val /= finalNorm;
                        MapPush(rowMap,colidx,val);
                    }
                    #endif
                }
                break;
            }
            case SPIN_DOWN:{
                if (cpcm(siteI, siteJ, pairRepI.second, pairRepI.first, pairRepIf.second, sign)){
                    #ifdef DISTRIBUTED_BASIS
                    MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                    #else
                    if (pt_Basis->search(pairRepIf, colidx)){
                        dataType val = factor * sign;
                        double finalNorm = pt_Basis->getNorm(colidx);
                        val /= finalNorm;
                        MapPush(rowMap,colidx,val);
                    }
                    #endif
                }
                break;
            }
        }
    };

    // void cmcp(SPIN spin, int siteI, int siteJ, dataType factor, pairIndex pairInd, int* initVec, int* initVecp, Basis* pt_Basis, MAP* rowMap){
    //     cpcm(spin, siteJ, siteI, factor, pairInd, initVec, initVecp, pt_Basis, rowMap);
    // };
};
class SpinOperator{
/*
    idx:0 --> spin:s
    idx++ ---> s--
*/
protected:
    Basis* pt_Basis;
    LATTICE_MODEL model;
    int spinDim;
    std::vector<double> szMat, spMat, smMat;
public: 
    SpinOperator(Basis* pt_Ba, LATTICE_MODEL mod=HEISENBERG);
    ~SpinOperator(){};
    
    double getSz(int siteI, VecI& initVec) const {return szMat.at(initVec.at(siteI));}
    void szsz(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        #ifdef DISTRIBUTED_BASIS
            dataType dval = factor * szMat.at(initVec[siteI]) * szMat.at(initVec[siteJ]);
            MapPush(rowMap,initInd,dval);
        #else
            ind_int colID;
            if (pt_Basis->search(initInd,colID)){
                dataType dval = factor * szMat.at(initVec[siteI]) * szMat.at(initVec[siteJ]);
                double finalNorm = pt_Basis->getNorm(colID);
                dval /= finalNorm;
                MapPush(rowMap,colID,dval);
            }
        #endif
    };
    
    void spsm(int siteI, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        #ifdef DISTRIBUTED_BASIS
            if (initVec[siteI] < (spinDim-1)){
                dataType val = factor * spMat[initVec[siteI]+1] * smMat[initVec[siteI]];
                MapPush(rowMap,initInd,val);
            }
        #else
            ind_int colID;
            if (pt_Basis->search(initInd, colID)){
                if (initVec[siteI] < (spinDim-1)){
                    dataType val = factor * spMat[initVec[siteI]+1] * smMat[initVec[siteI]];
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            }
        #endif
    };
    
    void smsp(int siteI, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        #ifdef DISTRIBUTED_BASIS
            if (initVec[siteI] > 0){
                dataType val = factor * smMat[initVec[siteI]-1] * spMat[initVec[siteI]];
                MapPush(rowMap,initInd,val);
            }
        #else
            ind_int colID;
            if (pt_Basis->search(initInd, colID)){
                if (initVec[siteI] > 0){
                    dataType val = factor * smMat[initVec[siteI]-1] * spMat[initVec[siteI]];
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            }
        #endif
    }
    
    void spsm(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        if (siteI==siteJ){
            spsm(siteI, factor, initInd, initVec, rowMap);
            return;
        }
        if (initVec[siteI] > 0 and initVec[siteJ] < (spinDim-1)){
            ind_int finalInd = initInd - pt_Basis->getmul(siteI) + pt_Basis->getmul(siteJ);
            #ifdef Distributde_Basis
                dataType val = factor * spMat[initVec[siteI]] * smMat[initVec[siteJ]];
                MapPush(rowMap,initInd,val);
            #else
                ind_int colID;
                if (pt_Basis->search(finalInd, colID)){
                    dataType val = factor * spMat[initVec[siteI]] * smMat[initVec[siteJ]];
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            #endif
        }
    };
    
    void smsp(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        if (siteI==siteJ){
            smsp(siteI, factor, initInd, initVec, rowMap);
            return;
        }
        spsm(siteJ, siteI, factor, initInd, initVec,rowMap);
    };
    /*
        Binary Reps For Spindim=2
    */
   double getSz(int siteI, ind_int repI) const {return szMat.at(1&(repI>>siteI));}
    void szsz(int siteI, int siteJ, dataType factor, ind_int repI, MAP* rowMap){
        #ifdef DISTRIBUTED_BASIS
            dataType dval = factor * szMat.at(1&(repI>>siteI)) * szMat.at(1&(repI>>siteJ));
            MapPush(rowMap,repI,dval);
        #else
            ind_int colID;
            if (pt_Basis->search(repI,colID)){
                dataType dval = factor * szMat.at(1&(repI>>siteI)) * szMat.at(1&(repI>>siteJ));
                double finalNorm = pt_Basis->getNorm(colID);
                dval /= finalNorm;
                MapPush(rowMap,colID,dval);
            }
        #endif
    };
    
    void spsm(int siteI, dataType factor, ind_int repI, MAP* rowMap){
        if (!bitTest(repI,siteI)){
            #ifdef DISTRIBUTED_BASIS
                dataType val = factor * spMat[1] * smMat[0];
                MapPush(rowMap,repI,val);
            #else
                ind_int colID;
                if (pt_Basis->search(repI, colID)){
                    dataType val = factor * spMat[1] * smMat[0];
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            #endif
        }
    };
    
    void smsp(int siteI, dataType factor, ind_int repI, MAP* rowMap){
        if (bitTest(repI,siteI)){
            #ifdef DISTRIBUTED_BASIS
                dataType val = factor * spMat[1] * smMat[0];
                MapPush(rowMap,repI,val);
            #else
                ind_int colID;
                if (pt_Basis->search(repI, colID)){
                    dataType val = factor * spMat[1] * smMat[0];
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            #endif
        }
    }
    
    void spsm(int siteI, int siteJ, dataType factor, ind_int repI, MAP* rowMap){
        ind_int colID;
        if (siteI==siteJ){
            spsm(siteI, factor, repI, rowMap);
            return;
        }
        if (bitTest(repI,siteI) && (!bitTest(repI,siteJ))){
            bitFlip(repI,siteI);
            bitFlip(repI,siteJ);
            dataType val = factor * spMat[1] * smMat[0];
            #ifdef DISTRIBUTED_BASIS
                MapPush(rowMap,repI,val);
            #else
                ind_int colID;
                if (pt_Basis->search(repI, colID)){
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
            #endif
        }
    };
    
    void smsp(int siteI, int siteJ, dataType factor, ind_int repI, MAP* rowMap){
        if (siteI==siteJ){
            smsp(siteI, factor, repI,rowMap);
            return;
        }
        spsm(siteJ, siteI, factor, repI,rowMap);
    };
};

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
        FermionOperator(pt_Ba),SparseMatrix<dataType>(pt_Ba, pt_Ba, pt_Ba->getSubDim(),spmNum_,dmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
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
    Hubbard& pushV(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++)V.at(pt_lattice->getOrbID(*it))=val;return *this;}
    Hubbard& pushU(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++)U.at(pt_lattice->getOrbID(*it))=val;return *this;}
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
public:
    Current(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1):FermionOperator(pt_Ba),SparseMatrix<cdouble>(pt_Ba, pt_Ba, pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
    ~Current(){};
    Current& pushLink(Link<cdouble> link, int matID){
        Links.push_back(link); Links[linkCount].setid(linkCount,matID); Links[linkCount].genLinkMaps(pt_lattice); 
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
        SpinOperator(pt_Ba,HEISENBERG),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
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
        SpinOperator(pt_Bi,HEISENBERG),SparseMatrix<cdouble>(pt_Bi,pt_Bf,pt_Bf_->getSubDim(),spmNum_){
            assert(pt_Bi->getPGIndex()==-1 and pt_Bf->getPGIndex()==-1);
            Ki = pt_Bi->getkIndex();
            Kf = pt_Bf->getkIndex();
            // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
            for (int i = 0; i < pt_lattice->getSiteNum(); i++) expFactor[i] = pt_lattice->expKR(Kf,i)/pt_lattice->expKR(Ki,i);
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
    double val = diagVal(occ,docc);
    SparseMatrix<T>::putDiag(val,rowID);
    // off diagonal part
    std::vector<ind_int> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); i++){
        pairIndex pairRepI = pt_Basis->getPairRepI(finalIndList[i]);
        bool isfminRep = pt_Basis->isfMin(pairRepI.first);
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
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
void Heisenberg<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    if(pt_Basis->getSiteDim()==2){
        int kIndex = pt_Basis->getkIndex();
        double initNorm, finalNorm;
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        initNorm = pt_Basis->getNorm(rowID);
        for (int i = 0; i < finalIndList.size(); i++){
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                int matID = (*linkit).getmatid();
                cdouble factor = factorList.at(i) * (*linkit).getVal();
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
            }
        }
    }
    else{
        int kIndex = pt_Basis->getkIndex();
        VecI initVec(pt_lattice->getOrbNum());
        double initNorm, finalNorm;
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        initNorm = pt_Basis->getNorm(rowID);
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
    double initNorm, finalNorm;
    // calculate <R1k|H*Pk|R2k>/norm1/norm2
    for (ind_int rowID = SparseMatrix<T>::startRow; rowID < SparseMatrix<T>::endRow; rowID++){
        rowMap.clear();
        std::vector<cdouble> factorList;
        std::vector<ind_int> finalIndList;
        pt_Basis->genSymm(rowID, finalIndList, factorList);
        initNorm = pt_Basis->getNorm(rowID);
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
    SpinOperator(pt_Ba,HEISENBERG),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_){
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
        double initNorm, finalNorm;
        initNorm = pt_Basis->getNorm(rowID);
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
        double initNorm, finalNorm;
        initNorm = pt_Basis->getNorm(rowID);
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
    double initNorm, finalNorm;
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        initNorm = pt_Basis->getNorm(rowID);
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
                ind_int repI = pt_Bi->getRepI(colID);
                cdouble dval = 0.0;
                for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                    dval += getSz(siteID,repI) * expFactor[siteID];
                }
                dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
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

#endif
