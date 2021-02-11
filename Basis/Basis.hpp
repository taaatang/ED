//
//  Basis.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Basis_hpp
#define Basis_hpp

#include "Global/globalPara.hpp"
#include "Geometry/Geometry.hpp"
#include "Utils/utils.hpp"

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
    Basis():model(LATTICE_MODEL::HUBBARD), kIndex(-1), PGRepIndex(-1){};
    // occList contains occupation number on each site dimension
    Basis(LATTICE_MODEL input_model, Geometry *pt_lat, VecI occList, int kInd=-1, int PGRepInd = -1, int siteD=2);
    ~Basis(){};

    LATTICE_MODEL getModel() const {return model;}
    int getkIndex() const {return kIndex;}
    int getPGIndex() const {return PGRepIndex;}

    // Rep integer access
    // rowidx->repI
    idx_t getRepI(idx_t idx) const {if(kIndex==-1 and model==LATTICE_MODEL::HUBBARD)return idx; return indexList.at(idx);}
    // repI->pairRepI
    pairIndex getPairRepI(idx_t repI) const {assert_msg(model!=LATTICE_MODEL::HEISENBERG,"PairRepI not defined for Heisenberg Model!");return std::make_pair(fIndexList.at(repI/sDim), sIndexList.at(repI%sDim));}
    // pairRepI->repI
    idx_t getRepI(pairIndex pairRepI) const {return fRepIdxHash.at(pairRepI.first)*sDim+sRepIdxHash.at(pairRepI.second);}
    
    idx_t getmul(int orbid) const {return mul.at(orbid);}

    void initMinMaxRep() const;

    //generate fIdex,sIdex list for Hubbard and tJ model
    void gendcmp();
    
    // construc k-subspace basis. need lattice to know symmetry operations
    void gen();
    void gen(int workerID, int workerNum);
    
    // construct k-subspace basis/norm from reps loaded from file
    void gen(std::string basisfile);
    void gen(std::string basisfile, std::string normfile);
    // distributed basis
    void gen(std::string basisfile, std::string normfile, int workerID, int workerNum);
    
    int getSiteDim() const {return siteDim;}
    int getOrbNum() const {return N;}
    idx_t getTotDim() const {return totDim;}
    idx_t getSubDim() const {return subDim;}
    idx_t getLocDim() const {return locDim;}
    
    double getNorm(idx_t rowid) const {
        if (kIndex==-1) return 1.0;
        #ifdef KEEP_BASIS_NORM
            return normList.at(rowid);
        #else
            return Norm(getRepI(rowid));
        #endif
    };
    
    // vec to Basis index
    idx_t vecToRep(VecI& v) const;
    pairIndex vecToRep(VecI& v, VecI& vp) const;
    
    // Basis index to Basis vec
    void repToVec(idx_t index, VecI& v) const;
    void repToVec(pairIndex pairInd, VecI& v, VecI& vp) const;
    void repToVec(idx_t index, VecI& v, VecI& vp) const;

    // for distributed basis. get the workerID where repI is possiblely stored
    bool getBid(idx_t repI, int &bid) const;

    // Binary search the position of index in indexList

    //search full hilbert space
    idx_t search(idx_t repI, const std::vector<idx_t> &indList) const;
    idx_t search(idx_t repI) const;
    idx_t search(pairIndex pairRepI) const;
    //search current symm sector subspace
    bool search(idx_t repI, idx_t &idx, const std::vector<idx_t> &indList) const;
    bool search(idx_t repI, idx_t &idx) const;
    bool search(pairIndex pairRepI, idx_t &idx) const;

    /*
        ************
        * Symmetry *
        * **********
    */
    // apply all translation operations to a Basis vec indexed by ind.
    // finalInd contains all resulting basis indexes.
    void genSymm(idx_t ind, std::vector<idx_t>& finalInd) const;
    void genSymm(idx_t ind, std::vector<idx_t>& finalInd, std::vector<cdouble>& factorList) const;
    
    // judge if repI is min repI. append symm operations tp symmList if symm(repI)==repI. not guanranteed to be MinRep since its norm might = 0
    bool isMin(idx_t repI, VecI& symmList);
    bool isfMin(idx_t frepI) const {return kIndex==-1?true:fMinRepSymmHash.find(frepI)!=fMinRepSymmHash.end();}
    // judge if a basis vec belongs to k-subspace rep vecs
    // rep: smallest index and norm!=0
    bool isMinRep(idx_t repI, double& norm) const;
    
    // calculate the norm of a genery repI
    double Norm(idx_t repI) const;
    // calculate the norm of a minimum repI in a symm cycle. only defined for Hubbard and tJ model.
    double minNorm(idx_t repI) const;
    
    // I/O
    void print(std::ostream& os = std::cout) const;
    void saveBasis(std::string basisfile, bool is_app=false);
    void saveBasis(std::string basisfile, std::string normfile, bool is_app=false);

private:
    /*
     Instantiate an object managing the Basis set
    */
    Geometry *pt_lattice;
    // option: 0->Hubbard, 1->t-J, 2->Heisenberg
    LATTICE_MODEL model;

    int kIndex; // default initial value kIndex=-1 ==> full hilbertspace
    int PGRepIndex; // default -1 -> do not use point group symm

    int siteDim{2}; // Hilbert Space dimension of a single Site
    idx_t totDim{0}; // Total Hilbert Space dimension
    idx_t subDim{0}; // Total subspace dimension of corresponding symm
    idx_t locDim{0}; // Local stored dimension
    idx_t fDim{0}, sDim{0};

    // Heisenberg:Sztot in unit of hbar/2
    int Sztot{0};
    // Hubbard/tJ: N sites for spin-up, Np sites for spin-down.
    int N{0};
    // Nocc[i] is the occupation number at site dimension i.
    VecI Nocc; 
    /*
        HUBBARD/tJ Model: fIndexList->spin-up, sIndexList->spin-down
        repI = fIdx * len(sIndexList) + sIdx.
        pairRepI = (fIndexList(fidx), sIndexList(sidx))
    */
    std::vector<idx_t> fIndexList, sIndexList;
    std::unordered_map<idx_t,idx_t> fRepIdxHash,sRepIdxHash;
    std::vector<idx_t> fMinRepList;
    std::unordered_map<idx_t,VecI> fMinRepSymmHash;

    /* 
        repI list in corresponding subspace in ascending order
        For Hubbard, if no symm is used, this is empty.
    */
    std::vector<idx_t> indexList;
    std::unordered_map<idx_t,idx_t> repHashTable;
    // repIStratList[workerID] is the smallest repI stored in workerID's RAM
    std::vector<idx_t> repIStartList; 
    // repIEndList[workerID] is the largest repI stored in workerID's RAM
    std::vector<idx_t> repIEndList;
    // norm for reps in symmetrized indexList
    std::vector<double> normList; 
    
    // used for siteDim=2. Binary Rep
    mutable idx_t fminRep{0},fmaxRep{0},sminRep{0},smaxRep{0};
    /*
        used only if siteDim>2
        Heisenberg model: vec: 0 for spin-up 1 for spin-down
        Hubbard/tJ:vec:1 for spin-up electron. vecp: 1 for spin-down electrons
    */
    mutable std::vector<int> vec, vecp; // mutable:can be modified in const function
    // Multipliers. Basis vec index = sum(i):vec(i)*mul(i)
    std::vector<idx_t> mul;
};

#endif // Basis_hpp
