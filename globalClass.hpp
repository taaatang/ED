//
//  globalClass.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef globalClass_hpp
#define globalClass_hpp

#include "globalPara.hpp"
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


/*
   ************************************
   *   Geometry Abstract Base Class   *
   ************************************
*/
struct Site{
    /*
     Describe a single site of a lattice
     id holds the index of this site among all sites
     coord is the coordinate in the lattice basis
     */
    int id; // 0 ~ N-1
    int orbital; // Used to identify orbital in multiband model
    BasisLat coord; // in unit of ax, ay
    BasisXY coordxy; // in x-y coordinates
};

struct kSite{
    int id;
    int orbital;
    BasisXY coord; // in unit of bx, by
    BasisXY coordxy; /// in x-y coordinates
};


class Geometry{
public:
    /*
     Instantiate an object describing the real space geometry of a lattice
     */
    int Nx, Ny, N; // Number of sites
    int N_enlg; // Number of sites in the enlarged lattice (Used for mapping Site interacting with boundary Site back to lattice).
    bool is_PBC; // Boundary condition
    std::vector<std::vector<int>> TransList;
    std::string name;
    
    BasisLat a1, a2; // Lattice basis
    BasisXY ax, ay; // The corresponding values of lattice basis in xy coordinates
    BasisXY b1, b2, kx, ky; // dual lattice basis
    
    Site* lattice_; // Pointer to a list of Site in the lattice
    Site* enlgLattice_; // Pointer to a list of Site in the enlarged lattice
    
    kSite* KLattice_; // Dual K space lattice
    
    BondMap *bondMaps_[MAX_BOND_NUM];
    std::string bondNames[MAX_BOND_NUM];
    int bondMapNum;
    
    Geometry();
    ~Geometry();
    
    std::string getName(); // Lattice name
    int getSiteNum();
    int* indToCoord(int i); // Site index -> coord in ax-ay coordinates
    double* indToCoordxy(int i); // Site index -> coord in x-y coordinates
    bool coordToInd(int* coord, int &ind); // ax-ay coord -> site index
    void genBonds(int linkNum, BasisLat *links, std::string bondName);
    // generate translation map for all r
    void genTransList();
    void print();
    
//    virtual void genKLattice() = 0; // Generate dual K space lattice;
};

class TriAngLattice: public Geometry{
public:
    // constructor for lattice with C6 Symm
    TriAngLattice(int numSites);
    // constructor for N1*N2 lattice
    TriAngLattice(int N1, int N2);
    ~TriAngLattice();
};

class SquareLattice: public Geometry{};

/*
    ***************
    * Basis Class *
    ***************
*/
class Basis{
public:
    /*
     Instantiate an object managing the Basis set
    */
    Geometry *pt_lattice;
    int kIndex; // default initial value kIndex=-1 ==> full hilbertspace
    int siteDim; // Hilbert Space dimension of a single Site
    ind_int totDim; // Total Hilbert Space dimension
    ind_int subDim; // Total subspace dimension
    int N; // Number of lattice sites
    
    std::vector<ind_int> indexList; // A list containing the index of basis vectors in an ascending order
    std::vector<double> normList;
    std::vector<ind_int> repIndList; // reps for a translation symmetry cycle (smallest index)
    
    int* initVec_;
    int* finalVec_;
    ind_int* finalInd;
    
    // Basis vec. Heisenberg model
    std::vector<int> vec;
    // Basis vec for spin-up and spin-down. Hubbard and t-J model
    std::vector<int> uvec, dvec;
    // Multipliers. Basis vec index = sum(i):vec(i)*mul(i)
    std::vector<ind_int> mul;
    
    // dimList contains occupation number on each site dimension
    Basis(int siteDim, int* dimList, Geometry *pt_lat);
    ~Basis();
    // generate full hilbertspace basis;
    void gen();
    // construc k-subspace basis. need lattice to know symmetry operations
    void gen(int kInd);
    // construct k-subspace basis/norm from reps loaded from file
    void gen(int kInd, std::string basisfile);
    void gen(int kInd, std::string basisfile, std::string normfile);
    // generate reps for translation symm
    void genRep();
    
    int getSiteDim();
    
    ind_int getTotDim();
    
    double getNorm(ind_int ind){
        if (kIndex==-1) return 1.0;
    #ifdef KEEP_BASIS_NORM
        return normList.at(ind);
    #else
        return kNorm(kIndex, indexList.at(ind), pt_lattice);
    #endif
    };
    
    // vec to Basis index
    ind_int vecToInd(int* v);
    
    // Basis index to Basis vec
    void indToVec(ind_int index, int* output);
    
    // Binary search the position of index in indexList
    ind_int search(ind_int index);
    bool search(ind_int index, ind_int &ind);
    
    // apply all translation operations to a Basis vec indexed by ind.
    // finalInd contains all resulting basis indexes.
    void genTranslation(Geometry* pt_lattice, ind_int ind, ind_int* finalInd);
    
    // judge if a basis vec belongs to k-subspace rep vecs
    // rep: smallest index and norm!=0
    bool isKRep(int kInd, ind_int initInd, Geometry* pt_lattice, double& norm);
    
    // calculate the norm of |rk>
    double kNorm(int kInd, ind_int initInd, Geometry* pt_lattice);
    
    // I/O
    void saveBasis(std::string basisfile);
    void saveBasis(std::string basisfile, std::string normfile);
};
    
 /* globalClass_hpp */

/*
    *******************
    * Operators Class *
    *******************
*/
class SpinOperator{
public:
    int spinDim;
    std::vector<double> szMat, spMat, smMat;
    SpinOperator();
    SpinOperator(int dim);
    ~SpinOperator();
    
    void szsz(int siteI, int siteJ, dataType factor, ind_int initInd, int* initVec, Basis* pt_Basis, MAP* rowMap){
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
    
    void spsm(int siteI, dataType factor, ind_int initInd, int* initVec, Basis* pt_Basis, MAP* rowMap){
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
    
    void smsp(int siteI, dataType factor, ind_int initInd, int* initVec, Basis* pt_Basis, MAP* rowMap){
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
    
    void spsm(int siteI, int siteJ, dataType factor, ind_int initInd, int* initVec, Basis* pt_Basis, MAP* rowMap){
        ind_int colID;
        if (siteI==siteJ){
            spsm(siteI, factor, initInd, initVec, pt_Basis, rowMap);
            return;
        }
        if (initVec[siteI] > 0 and initVec[siteJ] < (spinDim-1)){
            ind_int finalInd = initInd - pt_Basis->mul[siteI] + pt_Basis->mul[siteJ];
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
    
    void smsp(int siteI, int siteJ, dataType factor, ind_int initInd, int* initVec, Basis* pt_Basis, MAP* rowMap){
        if (siteI==siteJ){
            smsp(siteI, factor, initInd, initVec, pt_Basis, rowMap);
            return;
        }
        spsm(siteJ, siteI, factor, initInd, initVec, pt_Basis, rowMap);
    };
};

#ifdef SAXPY
    class HeisenbergOneHalf: public SpinOperator, public SparseMatrixSaxpy<dataType>{
#else
    class HeisenbergOneHalf: public SpinOperator, public SparseMatrix<dataType>{
#endif
    public:
    #ifdef SAXPY
        std::vector<ind_int> colListJ1, colListJ2, rowColInitListJ1, rowColInitListJ2;
        std::vector<dataType> valListJ1, valListJ2;
    #endif
        HeisenbergOneHalf(ind_int totDim, dataType J1, dataType J2);
        ~HeisenbergOneHalf();
        
        void genMat(Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum);
        void genMatOffDiag(Geometry *pt_lattice, Basis *pt_Basis, int couplingID, int polarNum);
        void genMatPol(Geometry *pt_lattice, Basis *pt_Basis, int couplingID, int polarID);
        void genDiagMat(Geometry *pt_lattice, Basis *pt_Basis, int couplingID, int polarNum);
        
        void genSubMat(int kIndex, Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum);
        void genSubMatMap(int kIndex, Geometry *pt_lattice, Basis *pt_Basis, int couplingNum, int polarNum);
        void setJ2(dataType J2);
        
        template<class T>
        void save(std::vector<T> *pt_vec, std::ofstream *pt_file, std::string filename);
        
        void saveAll(std::string dataDir);
        void HdiagPARPACK();
    };

class SzqOneHalf: public SpinOperator, public SparseMatrix<cdouble>{
public:
    SzqOneHalf(){};
    SzqOneHalf(ind_int totDim);
//    SzqOneHalf(ind_int rowDim, ind_int colDim);
    ~SzqOneHalf();
    
    void genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q);
    void genMat(Geometry* pt_lattice, Basis* pt_B1, Basis* pt_B2, BasisXY q);
};

class SSOneHalf: public SpinOperator, public SparseMatrix<dataType>{
public:
    SSOneHalf(){};
    SSOneHalf(ind_int totDim);
    ~SSOneHalf();
    
    void genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q, int initSiteID = 0);
    void genPairMat(Geometry* pt_lattice, Basis* pt_Basis, int siteJ, int siteI = 0);
    void genPairMat(int kInd, Geometry* pt_lattice, Basis* pt_Basis, int rIndex);
};

class comm_mpi{
    
};

class SparseEigenSolver{
    
};

class LanczosIterator{
    
};

class KrylovTimeEvolver{
    
};

class Spectra{
    
};

class Raman{
    						
};

#endif
