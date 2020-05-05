//
//  Basis.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Basis_hpp
#define Basis_hpp

#include "../Global/globalPara.hpp"
#include "../Geometry/Geometry.hpp"
#include "../Utils/utils.hpp"

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

/*
    ***************
    * Basis Class *
    ***************
*/
class Basis{
/*
    For Heisenberg model:
            repI is the binary represnetation of the spin configurations
    For Hubbard model,:
        repI=fidx*sDim+sidx. fIndexList[fidx] and sIndexList[sidx] are corresponding binary representations of
        spin-up and spin-down occupations.
        pairRepI = (fIndexList[fidx], sIndexList[sidx])
*/
public:  
    Basis():model(MODEL), kIndex(-1), PGRepIndex(-1){};
    // occList contains occupation number on each site dimension
    Basis(LATTICE_MODEL input_model, Geometry *pt_lat, VecI& occList, int kInd=-1, int PGRepInd = -1, int siteD=2);
    ~Basis(){};

    int getkIndex() const {return kIndex;}
    int getPGIndex() const {return PGRepIndex;}

    // Rep integer access
    // rowidx->repI
    ind_int getRepI(ind_int idx) const {if(kIndex==-1 and model==LATTICE_MODEL::HUBBARD)return idx; return indexList.at(idx);}
    // repI->pairRepI
    pairIndex getPairRepI(ind_int repI) const {return std::make_pair(fIndexList.at(repI/sDim), sIndexList.at(repI%sDim));}
    // pairRepI->repI
    ind_int getRepI(pairIndex pairRepI) const {return search(pairRepI);}
    // for Heisenberg. idx->repI
    void getRepI(ind_int idx, ind_int& repI) const {repI = indexList.at(idx);}
    // for Hubbard/tJ. repI->pairRepI
    void getRepI(ind_int repI, pairIndex& pairRepI) const {pairRepI.first = fIndexList.at(repI/sDim); pairRepI.second = sIndexList.at(repI%sDim);}
    
    ind_int getmul(int orbid) const {return mul.at(orbid);}

    void initMinMaxRep() const;

    //generate fIdex,sIdex list for Hubbard and tJ model
    void gendcmp();
    // generate full hilbertspace basis;
    void genFull();
    // construc k-subspace basis. need lattice to know symmetry operations
    void gen();
    
    // construct k-subspace basis/norm from reps loaded from file
    void gen(std::string basisfile);
    void gen(std::string basisfile, std::string normfile);
    // distributed basis
    void gen(std::string basisfile, std::string normfile, int workerID, int workerNum);
    
    int getSiteDim() const {return siteDim;}
    ind_int getTotDim() const {return totDim;}
    ind_int getSubDim() const {return subDim;}
    ind_int getLocDim() const {return locDim;}
    
    double getNorm(ind_int rowid) const {
        if (kIndex==-1) return 1.0;
        #ifdef KEEP_BASIS_NORM
            return normList.at(rowid);
        #else
            return Norm(kIndex, indexList.at(rowid), pt_lattice);
        #endif
    };
    
    // vec to Basis index
    ind_int vecToRep(VecI& v) const;
    pairIndex vecToRep(VecI& v, VecI& vp) const;
    
    // Basis index to Basis vec
    void repToVec(ind_int index, VecI& v) const;
    void repToVec(pairIndex pairInd, VecI& v, VecI& vp) const;
    void repToVec(ind_int index, VecI& v, VecI& vp) const;

    // for distributed basis. get the workerID where repI is possiblely stored
    bool getBid(ind_int repI, int &bid) const;

    // Binary search the position of index in indexList
    ind_int search(ind_int index, const std::vector<ind_int> &indList) const;
    ind_int search(ind_int index) const;
    ind_int search(pairIndex pairInd) const;//search full hilbert space
    bool search(ind_int index, ind_int &ind, const std::vector<ind_int> &indList) const;
    bool search(ind_int index, ind_int &ind) const;
    bool search(pairIndex pairInd, ind_int &ind) const;//search current symm sector subspace

    /*
        ************
        * Symmetry *
        * **********
    */
    // apply all translation operations to a Basis vec indexed by ind.
    // finalInd contains all resulting basis indexes.
    void genSymm(ind_int ind, std::vector<ind_int>& finalInd) const;
    void genSymm(ind_int ind, std::vector<ind_int>& finalInd, std::vector<cdouble>& factorList) const;
    
    // judge if a basis vec belongs to k-subspace rep vecs
    // rep: smallest index and norm!=0
    bool isMinRep(ind_int initInd, double& norm) const;
    
    // calculate the norm of |rk>
    double Norm(ind_int initInd) const;
    
    // I/O
    void saveBasis(std::string basisfile);
    void saveBasis(std::string basisfile, std::string normfile);
private:
    /*
     Instantiate an object managing the Basis set
    */
    Geometry *pt_lattice;
    // option: 0->Hubbard, 1->t-J, 2->Heisenberg
    LATTICE_MODEL model;

    int kIndex; // default initial value kIndex=-1 ==> full hilbertspace
    int PGRepIndex; // default -1 -> do not use point group symm

    int siteDim; // Hilbert Space dimension of a single Site
    ind_int totDim; // Total Hilbert Space dimension
    ind_int subDim; // Total subspace dimension of corresponding symm
    ind_int locDim; // Local stored dimension
    ind_int fDim, sDim;

    // Heisenberg:Sztot in unit of hbar/2
    int Sztot;
    // Hubbard/tJ: N sites for spin-up, Np sites for spin-down.
    int N;
    // Nocc[i] is the occupation number at site dimension i.
    VecI Nocc; 
    /*
        HUBBARD/tJ Model: fIndexList->spin-up, sIndexList->spin-down
        repI = fIdx * len(sIndexList) + sIdx.
        pairRepI = (fIndexList(fidx), sIndexList(sidx))
    */
    std::vector<ind_int> fIndexList, sIndexList;

    /* 
        repI list in corresponding subspace in ascending order
        For Hubbard, if no symm is used, this is empty.
    */
    std::vector<ind_int> indexList;
    // repIStratList[workerID] is the smallest repI stored in workerID's RAM
    std::vector<ind_int> repIStartList; 
    // repIEndList[workerID] is the largest repI stored in workerID's RAM
    std::vector<ind_int> repIEndList;
    // norm for reps in symmetrized indexList
    std::vector<double> normList; 
    
    // used for siteDim=2. Binary Rep
    mutable ind_int fminRep,fmaxRep,sminRep,smaxRep;
    /*
        used only if siteDim>2
        Heisenberg model: vec: 0 for spin-up 1 for spin-down
        Hubbard/tJ:vec:1 for spin-up electron. vecp: 1 for spin-down electrons
    */
    mutable std::vector<int> vec, vecp; // mutable:can be modified in const function
    // Multipliers. Basis vec index = sum(i):vec(i)*mul(i)
    std::vector<ind_int> mul;
};

#endif // Basis_hpp
