//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Operators_hpp
#define Operators_hpp

#include "globalPara.hpp"
#include "Basis.hpp"
#include "Geometry.hpp"
#include "SparseMatrix.hpp"
#include "utils.hpp"

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
        double finalNorm = pt_Basis->getNorm(rowID);
        factor /= finalNorm; 
        auto it = rowMap->find(rowID);
        if (it == rowMap->end()){
            (*rowMap)[rowID] = factor;
        }else{
            it->second += factor;
        }
    }

    bool cpcm(int siteI, int siteJ, ind_int initInd, const VecI& initVec, const VecI& initVecp, ind_int& finalInd, int& sign){
        if (siteI==siteJ){
            if (initVec[siteI]==1){
                finalInd = initInd;
                sign = 1;
                return true;
            }
        }else{
            switch(model){
                case HUBBARD:{
                    if (initVec[siteI]==0 && initVec[siteJ]==1){
                        finalInd = initInd + pt_Basis->getmul(siteI) - pt_Basis->getmul(siteJ);
                        int counter = 0;
                        if (siteI < siteJ) {
                            for (int i = siteI + 1; i < siteJ; i++) if (initVec[i]==1) counter++;
                        }else{
                            for (int i = siteJ + 1; i < siteI; i++) if (initVec[i]==1) counter++;
                        }
                        if (counter%2==0) sign = 1;
                        else sign = -1;
                        return true;
                    }
                    break;
                }

                // no double occupancy
                case t_J:{
                    if (initVec[siteI]==0 && initVecp[siteI]==0 && initVec[siteJ]==1){
                        finalInd = initInd + pt_Basis->getmul(siteI) - pt_Basis->getmul(siteJ);
                        int counter = 0;
                        if (siteI < siteJ) {
                            for (int i = siteI + 1; i < siteJ; i++) if (initVec[i]==1) counter++;
                        }else{
                            for (int i = siteJ + 1; i < siteI; i++) if (initVec[i]==1) counter++;
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
    void cpcm(SPIN spin, int siteI, int siteJ, dataType factor, ind_int repInd, const VecI& initVec, const VecI& initVecp, MAP* rowMap){
        pairIndex pairInd=pt_Basis->getpair(repInd);
        int sign;
        ind_int finalInd;
        ind_int colID;
        switch(spin){
            case SPIN_UP:{
                if (cpcm(siteI, siteJ, pairInd.first, initVec, initVecp, finalInd, sign)){
                    if (pt_Basis->search(std::make_pair(finalInd,pairInd.second), colID)){
                        dataType val = factor * sign;
                        double finalNorm = pt_Basis->getNorm(colID);
                        val /= finalNorm;
                        auto it = rowMap->find(colID);
                        if (it == rowMap->end()){
                            (*rowMap)[colID] = val;
                        }else{
                            it->second += val;
                        }
                    }
                }
                break;
            }
            case SPIN_DOWN:{
                if (cpcm(siteI, siteJ, pairInd.second, initVecp, initVec, finalInd, sign)){
                    if (pt_Basis->search(std::make_pair(pairInd.first,finalInd), colID)){
                        dataType val = factor * sign;
                        double finalNorm = pt_Basis->getNorm(colID);
                        val /= finalNorm;
                        auto it = rowMap->find(colID);
                        if (it == rowMap->end()){
                            (*rowMap)[colID] = val;
                        }else{
                            it->second += val;
                        }
                    }
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
protected:
    Basis* pt_Basis;
    LATTICE_MODEL model;
    int spinDim;
    std::vector<double> szMat, spMat, smMat;
public: 
    SpinOperator(Basis* pt_Ba, LATTICE_MODEL mod=HEISENBERG, int dim = 2);
    ~SpinOperator(){};
    
    double getSz(int siteI, VecI& initVec) const {return szMat.at(initVec.at(siteI));}
    void szsz(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        ind_int colID;
        if (pt_Basis->search(initInd,colID)){
            dataType dval = factor * szMat.at(initVec[siteI]) * szMat.at(initVec[siteJ]);
            double finalNorm = pt_Basis->getNorm(colID);
            auto it = rowMap->find(colID);
            if (it == rowMap->end()){
                (*rowMap)[colID] = dval/finalNorm;
            }else{
                it->second += dval/finalNorm;
            }
        }
    };
    
    void spsm(int siteI, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        ind_int colID;
        if (pt_Basis->search(initInd, colID)){
            if (initVec[siteI] < (spinDim-1)){
                dataType val = factor * spMat[initVec[siteI]+1] * smMat[initVec[siteI]];
                double finalNorm = pt_Basis->getNorm(colID);
                auto it = rowMap->find(colID);
                if (it == rowMap->end()){
                    (*rowMap)[colID] = val/finalNorm;
                }else{
                    it->second += val/finalNorm;
                }
            }
        }
    };
    
    void smsp(int siteI, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        ind_int colID;
        if (pt_Basis->search(initInd, colID)){
            if (initVec[siteI] > 0){
                dataType val = factor * smMat[initVec[siteI]-1] * spMat[initVec[siteI]];
                double finalNorm = pt_Basis->getNorm(colID);
                auto it = rowMap->find(colID);
                if (it == rowMap->end()){
                    (*rowMap)[colID] = val/finalNorm;
                }else{
                    it->second += val/finalNorm;
                }
            }
        }
    }
    
    void spsm(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        ind_int colID;
        if (siteI==siteJ){
            spsm(siteI, factor, initInd, initVec, rowMap);
            return;
        }
        if (initVec[siteI] > 0 and initVec[siteJ] < (spinDim-1)){
            ind_int finalInd = initInd - pt_Basis->getmul(siteI) + pt_Basis->getmul(siteJ);
            ind_int colID;
            if (pt_Basis->search(finalInd, colID)){
                dataType val = factor * spMat[initVec[siteI]] * smMat[initVec[siteJ]];
                double finalNorm = pt_Basis->getNorm(colID);
                auto it = rowMap->find(colID);
                if (it == rowMap->end()){
                    (*rowMap)[colID] = val/finalNorm;
                }else{
                    it->second += val/finalNorm;
                }
            }
        }
    };
    
    void smsp(int siteI, int siteJ, dataType factor, ind_int initInd, VecI& initVec, MAP* rowMap){
        if (siteI==siteJ){
            smsp(siteI, factor, initInd, initVec,rowMap);
            return;
        }
        spsm(siteJ, siteI, factor, initInd, initVec,rowMap);
    };
};

template <class T>
class Heisenberg: public SpinOperator, public SparseMatrix<T>{
private:
    // constant links
    int linkCount;
    int spmCount;
    std::vector<Link<double>*> Links;
    // non-constant links
    std::vector<Link<double>*> NCLinks;
    Geometry *pt_lattice;
public:
    Heisenberg(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2):\
        SpinOperator(pt_Ba,HEISENBERG,spindim),SparseMatrix<T>(pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
    ~Heisenberg(){};

    Heisenberg& pushLink(Link<double>& link){
        if(link.isConst()) {Links.push_back(&link); link.setid(linkCount,0); linkCount++; SparseMatrix<T>::parameters.at(0) = 1.0; }
        else {NCLinks.push_back(&link); spmCount++; link.setid(linkCount,spmCount); linkCount++; SparseMatrix<T>::parameters.at(spmCount) = 1.0;}
        link.genLinkMaps(pt_lattice); 
        return *this;
    }

    void setVal(Link<double>& link, double val){link.setVal(val); parameters.at(link.getmatid()) = val;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    void genMat();
};

template <class T>
class Hubbard: public FermionOperator, public SparseMatrix<dataType>{
private:
    // constant links
    int linkCount;
    int spmCount;
    std::vector<Link<double>*> Links;
    // non-constant links
    std::vector<Link<double>*> NCLinks;
    // V:charge transfer energy on different orbitals. U:onsite repulsion
    VecD V, U;
    Geometry *pt_lattice;
public:
    Hubbard(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2):V(pt_lat->getUnitOrbNum()), U(pt_lat->getUnitOrbNum()),\
        FermionOperator(pt_Ba),SparseMatrix<dataType>(pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){}
    ~Hubbard(){};

    Hubbard& pushLink(Link<double>& link){
        if(link.isConst()) {Links.push_back(&link); link.setid(linkCount,0); linkCount++; parameters.at(0) = 1.0; }
        else {NCLinks.push_back(&link); spmCount++; link.setid(linkCount,spmCount); linkCount++; parameters.at(spmCount) = 1.0;}
        link.genLinkMaps(pt_lattice); 
        return *this;
    }
    Hubbard& pushV(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++)V.at(pt_lattice->getOrbID(*it))=val;return *this;}
    Hubbard& pushU(std::vector<ORBITAL> orbList, double val){for(auto it=orbList.begin();it!=orbList.end();it++)U.at(pt_lattice->getOrbID(*it))=val;return *this;}
    void printV() const {std::cout<<"V:"<<std::endl;for(auto it=V.begin();it!=V.end();it++)std::cout<<*it<<", ";std::cout<<std::endl;}
    void printU() const {std::cout<<"U:"<<std::endl;for(auto it=U.begin();it!=U.end();it++)std::cout<<*it<<", ";std::cout<<std::endl;}
    double diagVal(VecI occ, VecI docc) const {double val=0.0; for(int i=0;i<pt_lattice->getUnitOrbNum();i++)val += occ.at(i)*V.at(i)+docc.at(i)*U.at(i);return val;}
    void setVal(Link<double>& link, double val){link.setVal(val); parameters.at(link.getmatid()) = val;}
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
        SpinOperator(pt_Bi,HEISENBERG,spindim),SparseMatrix<cdouble>(pt_Bf_->getSubDim(),spmNum_){
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
    int r;
    std::vector<VecI> siteJList;
public:
    SSOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2);
    ~SSOp(){}
    void setr(int r_){assert(r_<pt_lattice->getSiteNum());r=r_;}
    void row(ind_int rowID, std::vector<MAP>& rowMaps);
    // S(i)*S(i+r)
    void genPairMat(int rIndex);
    void project(double s, std::vector<T>& vec);
};

/*
    *********************
    * Hamiltonian Class *
    *********************
*/
template <class T>
void Heisenberg<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    int kIndex = pt_Basis->getkIndex();
    VecI initVec(pt_lattice->getOrbNum());
    double initNorm, finalNorm;
    std::vector<ind_int> finalIndList;
    pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
    initNorm = pt_Basis->getNorm(rowID);
    for (int i = 0; i < finalIndList.size(); i++){
        pt_Basis->indToVec(finalIndList[i], initVec);
        cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            cdouble factor1 = factor * (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, factor1, finalIndList[i], initVec, &rowMaps[0]);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMaps[0]);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMaps[0]);
            }
        }
    }

    for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
        int matID = (*linkit)->getmatid();
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
            factor *= (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
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
        std::vector<ind_int> finalIndList;
        pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
        initNorm = pt_Basis->getNorm(rowID);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor1 = factor * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor1, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMap);
                }
            }
        }
        SparseMatrix<T>::pushRow(&rowMap);

        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->indToVec(finalIndList[i], initVec);
                cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
                factor *= (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
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

template <class T>
void Hubbard<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    VecI initVec(pt_lattice->getOrbNum()), initVecp(pt_lattice->getOrbNum());
    /*
        *****************
        * Constant Part *
        * ***************
    */
    // diagonal part. occupancy and double-occ
    VecI occ, docc;
    ind_int initInd = pt_Basis->getRepI(rowID);
    pt_Basis->indToVec(initInd, initVec, initVecp);
    pt_lattice->orbOCC(initVec, initVecp, occ, docc);
    double val = diagVal(occ,docc);
    cdouble diag_val = 0.0;
    // off diagonal part
    std::vector<ind_int> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genTranslation(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); i++){
        pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
        if(finalIndList[i]==initInd)diag_val += val*factorList[i];
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            cdouble factor = factorList.at(i) * (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteI = (*bondit).at(0);
                int siteJ = (*bondit).at(1);
                // cp.siteI * cm.siteJ
                cpcm(SPIN::SPIN_UP, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMaps[0]);
                cpcm(SPIN::SPIN_UP, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMaps[0]);
                cpcm(SPIN::SPIN_DOWN, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMaps[0]);
                cpcm(SPIN::SPIN_DOWN, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMaps[0]);   
            }
        }
    }
    diag(rowID,diag_val,&rowMaps[0]);

    /*
        ***********************
        * Time Dependent Part *
        * *********************
    */
    for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
        int matID = (*linkit)->getmatid();
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
            cdouble factor = factorList.at(i) * (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                // cm.siteID * cp.siteIDp
                cpcm(SPIN::SPIN_UP, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMaps[matID]);
                cpcm(SPIN::SPIN_UP, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMaps[matID]);
                cpcm(SPIN::SPIN_DOWN, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap[matID]);
                cpcm(SPIN::SPIN_DOWN, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMaps[matID]); 
            }
        }
    }
}

template <class T>
void Hubbard<T>::genMat(){
    int kIndex = pt_Basis->getkIndex();
    clear();
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < spmNum; i++) pushRow(&rowMap,i);
    VecI initVec(pt_lattice->getOrbNum()), initVecp(pt_lattice->getOrbNum());
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        /*
            *****************
            * Constant Part *
            * ***************
        */
        rowMap.clear();
        // diagonal part. occupancy and double-occ
        VecI occ, docc;
        ind_int initInd = pt_Basis->getRepI(rowID);
        pt_Basis->indToVec(initInd, initVec, initVecp);
        pt_lattice->orbOCC(initVec, initVecp, occ, docc);
        double val = diagVal(occ,docc);
        cdouble diag_val = 0.0;
        // off diagonal part
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genTranslation(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
            if(finalIndList[i]==initInd)diag_val += val*factorList[i];
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor = factorList.at(i) * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
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
        pushRow(&rowMap);

        /*
            ***********************
            * Time Dependent Part *
            * *********************
        */
        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
                cdouble factor = factorList.at(i) * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // cm.siteID * cp.siteIDp
                    cpcm(SPIN::SPIN_UP, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_UP, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap); 
                }
            }
            pushRow(&rowMap, (*linkit)->getmatid());
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
    SpinOperator(pt_Ba,HEISENBERG,spindim),SparseMatrix<T>(pt_Ba->getSubDim(),spmNum_){
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
    int kIndex = pt_Basis->getkIndex();
    VecI initVec(pt_lattice->getOrbNum());
    double initNorm, finalNorm;
    initNorm = pt_Basis->getNorm(rowID);
    std::vector<ind_int> finalIndList;
    pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
    for (int i = 0; i < finalIndList.size(); i++){
        pt_Basis->indToVec(finalIndList[i], initVec);
        cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(pt_Basis->getkIndex(),i)/pt_lattice->getSiteNum()/initNorm;
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
template <class T>
void SSOp<T>::project(double s, std::vector<T>& vec){
    assert(r==-1);
    std::vector<T> tmp(vec.size());
    double smin = double(pt_lattice.getOrbNum()%2)/2.0;
    double smax = double(pt_lattice.getOrbNum())/2.0+0.1:
    for(double si = smin; si<smax; si += 1.0){
        if (std::abs(si-s)>INFINITESIMAL){
            MxV(vec, tmp.data());
            double stot = si * (si+1);
            #pragma omp parallel for
            for(ind_int i =0; i < vec.size();i++) vec[i] = tmp[i]-stot*vec[i];
        }
    }
}
// generate matrix in subsapce labeled by kIndex for sum.r:Sr*Sr+dr, dr is labeled by rIndex
template <class T>
void SSOp<T>::genPairMat(int rIndex){
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
        pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(pt_Basis->getkIndex(),i)/pt_lattice->getSiteNum()/initNorm;
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
}

template <class T>
void SzkOp<T>::row(ind_int rowID, std::vector<MAP>& rowMaps){
    ind_int colID;
    if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
        VecI initVec(pt_lattice->getOrbNum());
        cdouble dval = 0.0;
        pt_Bi->indToVec(pt_Bi->getRepI(colID), initVec);
        for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
            dval += getSz(siteID,initVec) * expFactor[siteID];
        }
        dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
        rowMaps[0][colID] = dval;
    }
}
template <class T>
void SzkOp<T>::genMat(){
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
                    pt_Bi->indToVec(pt_Bi->getRepI(colID), initVec);
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
                    pt_Bf->indToVec(pt_Bf->getRepI(rowID), initVec);
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
}

#endif
