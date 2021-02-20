//
//  Operators.hpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef Operators_hpp
#define Operators_hpp

#include <iostream>
#include <cmath>
#include <algorithm>

#include "Operator/OperatorsBase.hpp"
#include "Utils/utils.hpp"
#include "Pulse/pulse.hpp"
#include "Solver/PARPACKSolver.hpp"

/***************
 * HAMILTONIAN *
 ***************/

template <typename T> 
class HamiltonianBase : public OperatorBase<T> {
public:
    HamiltonianBase( ) { }
    HamiltonianBase(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_ = 1, int dmNum_ = 0);
    ~HamiltonianBase( ) { }

    // Add onsite energy V
    void pushV(std::vector<ORBITAL> orbList, double val);

    // Add onsite Coulomb interaction U
    void pushU(std::vector<ORBITAL> orbList, double val);

    void printV(std::ostream& os) const;

    void printU(std::ostream& os) const;

    // Calculate the sum of onsite energy and Coulomb interaction
    double diagVal(const VecI& occ, const VecI& docc) const;
    
    virtual void setPeierls(Pulse* pulse = nullptr) { }

    virtual void printPeierls(std::ostream& os = std::cout) { }

    virtual void row(idx_t rowidx, std::vector<MAP<T>>& rowMaps) = 0;

    PARPACKSolver<T> solver;

private:
    VecD V, U;
};

template <LATTICE_MODEL M, typename T>
class Hamiltonian: public HamiltonianBase<T> {
public:
    Hamiltonian(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_ = 1, int dmNum_ = 0);
    ~Hamiltonian( ) { }

    void setPeierls(Pulse* pulse = nullptr);

    void printPeierls(std::ostream& os = std::cout);

    bool next( );

    void row(idx_t rowidx, std::vector<MAP<T>>& rowMaps);

private:
    Pulse* pulse{nullptr};
    VecD PeierlsOverlap;
};



class Current: public OperatorBase<cdouble> {
public:
    Current(Geometry *latt, Basis *Ba);
    ~Current( ) { }

    void setDirection(std::string plz);
    void row(idx_t rowID, std::vector<MAP<cdouble>>& rowMaps);

private:
    void setDirection(VecD d);

    VecD direction;
    std::vector<Link<cdouble>> myHoppingT;
    std::string plz;
};

class Nocc: public OperatorBase<double> {
public:
    Nocc(Geometry *latt, Basis *Ba);
    ~Nocc(){}
    
    template <typename T>
    friend void save(T *d_pt, int size, std::string filename, bool is_app);
    template <typename T>
    friend void save(T *d_pt, idx_t size, std::string filename, bool is_app);

    void row(idx_t rowID, std::vector<MAP<double>>& rowMaps) { };
    inline void row(idx_t rowID);
    void construct( );

    double count(ORBITAL orbital, dataType* vec);
    double count(int orbid, dataType* vec);
    void count(dataType* vec);

    void clear( );
    void save(std::string dir);

private:
    std::vector<std::vector<double>> records;
};

template <class T>
class CkOp: public OperatorBase<T> {
public:
    CkOp(Geometry *latt, Basis *Bi, Basis *Bf);
    ~CkOp( ) { }

    void set(LADDER pm, SPIN spin, Orbital orb);
    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);

private:
    LADDER pm;
    SPIN spin;
    Orbital orb;
    int Ki, Kf;
    std::vector<cdouble> expFactor;
    VecI posList;
};

template<class T>
class RamanOp: public OperatorBase<T> {
public:
    RamanOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int dmNum = 0);
    ~RamanOp( ) { }

    void setplz(VecD pIn_, VecD pOut_);
    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);

private:
    // polarization of in/out photons
    VecD pIn, pOut;
    std::vector<VecD> RamanWeight;
    bool NoRW;
};

// correlator:Si*Si+r
// can also be used to construct total spin operator: S^2
template <class T>
class SSOp: public OperatorBase<T> {
public:
    SSOp(Geometry *latt, Basis* Bi, int spmNum = 1, int dmNum = 0, int spindim = 2);
    ~SSOp( ) { }

    void setr(int r_);

    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);
 
    void project(double s, T* vec);

private:
    // if r = -1, sum.r sum.i Si*Si+r. total S^2
    int r; 
    std::vector<VecI> siteJList;
};

template <class T>
class SzkOp: public OperatorBase<T> {
public:
    SzkOp(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum = 1, int dmNum = 0, int spindim = 2);
    ~SzkOp( ) { }

    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);

private:
    // initial k index Ki and final k index Kf
    int Ki, Kf;
    std::vector<cdouble> expFactor;
};

/*
    ****************************
    * Operator Implementations *
    ****************************
*/

template<typename T>
HamiltonianBase<T>::HamiltonianBase(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_, int dmNum_):\
OperatorBase<T>(latt, Bi, Bf, spmNum_, dmNum_), V(latt->getUnitOrbNum(),0.0), U(latt->getUnitOrbNum(),0.0) {
    solver = PARPACKSolver<T>(this, 1);
}

template <typename T>
void HamiltonianBase<T>::pushV(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbID(orb);
        for (auto& id : ids) {
            V.at(id) = val;
        }
    }
}

template <typename T>
void HamiltonianBase<T>::pushU(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbID(orb);
        for (auto& id : ids) {
            U.at(id) = val;
        }
    }
}

template <typename T>
void HamiltonianBase<T>::printV(std::ostream& os) const {
    os<<"V: ";
    for (auto val : V) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <typename T>
void HamiltonianBase<T>::printU(std::ostream& os) const {
    os<<"U: ";
    for (auto val : U) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <typename T>
double HamiltonianBase<T>::diagVal(const VecI& occ, const VecI& docc) const {
    double val{0.0};
    for (int i = 0; i < this->latt->getUnitOrbNum(); ++i) {
        val += occ.at(i) * V.at(i) + docc.at(i) * U.at(i);
    }
    return val;
}

template <LATTICE_MODEL M, typename T>
Hamiltonian<M, T>::Hamiltonian(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum, int dmNum):HamiltonianBase<T>(latt, Bi, Bf, spmNum, dmNum) {
    
}

template <LATTICE_MODEL M, typename T>
void Hamiltonian<M, T>::setPeierls(Pulse* pulse) {
    if constexpr (M == LATTICE_MODEL::HUBBARD or M == LATTICE_MODEL::tJ) {
        this->pulse = pulse;
        if (pulse) {
            auto epol = pulse->getPol();
            normalize(epol);
            VecD overlap;
            for(const auto& link : this->hoppingT) {
                assert_msg(link.getLinkVecNum()==1, "In setPeierls, each hopping link should only have one dr!");
                auto dr = link.getvec(0);
                dr = this->latt->RtoRxy(dr);
                auto val = vdotv(epol,dr);
                overlap.push_back(val);
                if (val!=0.0){
                    PeierlsOverlap.push_back(val);
                }
            }
            std::sort(PeierlsOverlap.begin(), PeierlsOverlap.end());
            auto last = std::unique(PeierlsOverlap.begin(), PeierlsOverlap.end());
            PeierlsOverlap.resize(last - PeierlsOverlap.begin());
            for (int idx = 0; idx < (int)overlap.size(); ++idx) {
                auto low = std::lower_bound(PeierlsOverlap.begin(), PeierlsOverlap.end(), overlap[idx]);
                if (low != PeierlsOverlap.end() and overlap[idx] == *low) {
                    int matid = 1 + 2 * (low - PeierlsOverlap.begin());
                    auto& link = this->hoppingT.at(idx);
                    link.setmatid(matid);
                    link.setConst(false);
                    link.setOrdered(true);
                }
            }
            this->setSpmNum(1 + 2 * PeierlsOverlap.size());
            PeierlsOverlap.insert(PeierlsOverlap.begin(), 0.0);
        }
    } else {
        std::cout<<"Peierls substitution is not defined for "<<M<<"\n";
        return;
    }
}

template<LATTICE_MODEL M, typename T>
void Hamiltonian<M, T>::printPeierls(std::ostream& os) {
    if constexpr (M == LATTICE_MODEL::HUBBARD or M == LATTICE_MODEL::tJ) {
        for (int id = 0; id < (int)PeierlsOverlap.size(); ++id) {
            os<<"A*dr = "<<PeierlsOverlap[id]<<"\n";
            auto matid = (id == 0) ? id : 2 * id -1;
            for (const auto& link : this->hoppingT) {
                if (link.getmatid() == matid) {
                    link.print();
                }
            }
            os<<"\n";
            // if (id > 0) {
            //     os<<"A*dr = "<<-PeierlsOverlap[id]<<"\n";
            //     ++matid;
            //     for (const auto& link : this->hoppingT) {
            //         if (link.getmatid() == matid) {
            //             link.print();
            //         }
            //     } 
            //     os<<"\n";
            // }
        }
    }
}

template<LATTICE_MODEL M, typename T>
bool Hamiltonian<M, T>::next( ) {
    if (pulse) {
        auto A = pulse->getAa();
        for (int i = 1; i < (int)PeierlsOverlap.size(); ++i) {
            auto factor = std::exp(CPLX_I * A * PeierlsOverlap[i]);
            auto idx = 2 * i - 1;
            this->setVal(idx, factor);
            this->setVal(idx + 1, std::conj(factor));
        }
        return pulse->next();
    } else {
        return false;
    }
}

template <LATTICE_MODEL M, typename T>
void Hamiltonian<M, T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if constexpr (M == LATTICE_MODEL::HUBBARD) {
        // diagonal part. occupancy and double-occ
        VecI occ, docc;
        idx_t repI = this->Bf->getRepI(rowID);
        pairIdx_t pairRepI = this->Bf->getPairRepI(repI);
        this->latt->orbOCC(pairRepI, occ, docc);
        T val = this->diagVal(occ,docc);
        for (const auto& link : this->interBandU){
            // auto bonds = link.bond();
            for (const auto& bond : link.bond()){
                auto siteI = bond.at(0);
                auto siteJ = bond.at(1);
                assert(siteI!=siteJ);
                val += link.getVal()*(double)(bitTest(pairRepI.first,siteI)+bitTest(pairRepI.second,siteI))*(double)(bitTest(pairRepI.first,siteJ)+bitTest(pairRepI.second,siteJ));
            }
        }
        SparseMatrix<T>::putDiag(val,rowID);
        // off diagonal part
        std::vector<pairIdx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->Bf->genSymm(rowID, finalIndList, factorList);
        int i = 0;
        for (auto pairRepI : finalIndList) {
            bool isfminRep = this->Bf->isfMin(pairRepI.first);
            for (const auto& link:this->hoppingT) {
                int matID = link.getmatid();
                int matIDp = matID; 
                if (link.isOrdered()) matIDp++;
                cdouble factor = factorList.at(i) * link.getVal();
                for (const auto& bond : link.bond()) {
                    int siteI = bond.at(0);
                    int siteJ = bond.at(1);
                    cdouble phase = this->latt->twistPhase(siteI,siteJ);
                    cdouble phase_c = std::conj(phase);
                    // cp.siteI * cm.siteJ
                    this->cpcm(SPIN::UP, siteI, siteJ, phase*factor, pairRepI, &rowMaps[matID]);
                    this->cpcm(SPIN::UP, siteJ, siteI, phase_c*factor, pairRepI, &rowMaps[matIDp]);
                    if(isfminRep){
                        this->cpcm(SPIN::DOWN, siteI, siteJ, phase*factor, pairRepI, &rowMaps[matID]);
                        this->cpcm(SPIN::DOWN, siteJ, siteI, phase_c*factor, pairRepI, &rowMaps[matIDp]);   
                    }
                }
            }

            for (const auto& link:this->exchangeJ) {
                int matID = link.getmatid();
                cdouble factor = factorList.at(i) * link.getVal();
                for (const auto& bond : link.bond()){
                    int siteI = bond.at(0);
                    int siteJ = bond.at(1);
                    this->exchange(siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
                }
            }

            for (const auto& link:this->pairHoppingJ) {
                int matID = link.getmatid();
                cdouble factor = factorList.at(i) * link.getVal();
                for (const auto& bond : link.bond()){
                    int siteI = bond.at(0);
                    int siteJ = bond.at(1);
                    this->pairHopping(siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
                }
            }
            ++i;
        }
    } else if constexpr(M == LATTICE_MODEL::tJ) {
        // off diagonal part
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->Bf->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < (int)finalIndList.size(); ++i){
            pairIdx_t pairRepI = this->Bf->getPairRepI(finalIndList[i]);
            bool isfminRep = this->Bf->isfMin(pairRepI.first);
            for (const auto& link:this->hoppingT) {
                int matID = link.getmatid();
                int matIDp = matID; 
                if (link.isOrdered()) ++matIDp;
                cdouble factor = factorList.at(i) * link.getVal();
                for (const auto& bond : link.bond()) {
                    int siteI = bond.at(0);
                    int siteJ = bond.at(1);
                    cdouble phase = this->latt->twistPhase(siteI,siteJ);
                    cdouble phase_c = std::conj(phase);
                    // cp.siteI * cm.siteJ
                    this->cpcm(SPIN::UP, siteI, siteJ, phase*factor, pairRepI, &rowMaps[matID]);
                    this->cpcm(SPIN::UP, siteJ, siteI, phase_c*factor, pairRepI, &rowMaps[matIDp]);
                    if(isfminRep){
                        this->cpcm(SPIN::DOWN, siteI, siteJ, phase*factor, pairRepI, &rowMaps[matID]);
                        this->cpcm(SPIN::DOWN, siteJ, siteI, phase_c*factor, pairRepI, &rowMaps[matIDp]);   
                    }
                }
            } 

            for (const auto& link:this->superExchangeJ) {
                int matID = link.getmatid();
                cdouble factor = factorList.at(i) * link.getVal();
                for (const auto& bond : link.bond()) {
                    int siteI = bond.at(0);
                    int siteJ = bond.at(1);
                    // szi * szj - 1/4*ni*nj 
                    this->szsznn(siteI, siteJ, factor, finalIndList[i], &rowMaps[matID]);
                    // 1/2 * sm.siteID * sp.siteIDp
                    this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[matID]);
                    // 1/2 * sp.siteID * sm.siteIDp
                    this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[matID]); 
                }
            } 
        }
        // diag(rowID,val,&rowMaps[0]);
    } else if constexpr(M == LATTICE_MODEL::HEISENBERG) {
        if(this->Bf->getSiteDim()==2){
            std::vector<idx_t> finalIndList;
            std::vector<cdouble> factorList;
            this->Bf->genSymm(rowID, finalIndList, factorList);
            for (int i = 0; i < (int)finalIndList.size(); ++i){
                for (auto linkit = this->Links.begin(); linkit != this->Links.end(); linkit++){
                    int matID = (*linkit).getmatid();
                    cdouble factor = factorList.at(i) * (*linkit).getVal();
                    switch((*linkit).getLinkType()){
                        case LINK_TYPE::SUPER_EXCHANGE_J:{
                            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                                int siteID = (*bondit).at(0);
                                int siteIDp = (*bondit).at(1);
                                // sz.siteID * sz.siteIDp
                                this->szsz(siteID, siteIDp, factor, finalIndList[i], &rowMaps[matID]);
                                // 1/2 * sm.siteID * sp.siteIDp
                                this->spsm(siteID, siteIDp, factor/2.0, finalIndList[i], &rowMaps[matID]);
                                // 1/2 * sp.siteID * sm.siteIDp
                                this->smsp(siteID, siteIDp, factor/2.0, finalIndList[i], &rowMaps[matID]);
                            }
                            break;
                        }
                        case LINK_TYPE::CHIRAL_K:{
                            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                                int siteI = (*bondit).at(0);
                                int siteJ = (*bondit).at(1);
                                int siteK = (*bondit).at(2);
                                this->chiral(siteI, siteJ, siteK, factor, finalIndList[i], &rowMaps[matID]);
                            }
                            break;
                        }
                        default: break;
                    }
                    
                }
            }
        } else {
            VecI initVec(this->latt->getOrbNum());
            std::vector<idx_t> finalIndList;
            std::vector<cdouble> factorList;
            this->Bf->genSymm(rowID, finalIndList, factorList);
            for (int i = 0; i < (int)finalIndList.size(); ++i){
                this->Bf->repToVec(finalIndList[i], initVec);
                for (auto linkit = this->Links.begin(); linkit != this->Links.end(); linkit++){
                    int matID = (*linkit).getmatid();
                    cdouble factor = factorList.at(i) * (*linkit).getVal();
                    for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                        int siteID = (*bondit).at(0);
                        int siteIDp = (*bondit).at(1);
                        // sz.siteID * sz.siteIDp
                        this->szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMaps[matID]);
                        // 1/2 * sm.siteID * sp.siteIDp
                        this->spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
                        // 1/2 * sp.siteID * sm.siteIDp
                        this->smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
                    }
                }
            }
        }  
    } else {
        std::cout<<" Input Lattice Model is not defined!\n";
        exit(1);
    }
}


template <class T>
CkOp<T>::CkOp(Geometry *latt, Basis *Bi, Basis *Bf):\
OperatorBase<T>(latt, Bi, Bf), expFactor(latt->getSiteNum()) {
    assert(Bi->getPGIndex()==-1 and Bf->getPGIndex()==-1);
    Ki = Bi->getkIndex();
    Kf = Bf->getkIndex();
    // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
    int N = latt->getSiteNum();
    for (int i = 0; i < N; ++i) {
        expFactor[i] = latt->expKR(Ki,i)/latt->expKR(Kf,i)/std::sqrt(N);
    }
}

template <class T>
void CkOp<T>::set(LADDER pm, SPIN spin, Orbital orb) {
    this->pm = pm;
    this->spin = spin;
    this->orb = orb;
    auto nis = this->Bi->getOcc();
    auto nfs = this->Bf->getOcc();
    auto ni = (spin == SPIN::UP) ? nis.at(0) : nis.at(1);
    auto nf = (spin == SPIN::UP) ? nfs.at(0) : nfs.at(1);
    auto dn = (pm == LADDER::PLUS) ? 1 : -1;
    assert_msg(ni + dn == nf, "CkOp init and final hilbert space particle numbers mismatch!");
    auto Norb = this->latt->getOrbNum();
    posList.clear();
    for (int i = 0; i < Norb; ++i) {
        if (this->latt->is_Orbital(i, orb.orb, orb.orbid)) posList.push_back(i);
    }
    assert((int)posList.size() == this->latt->getSiteNum());
}

template <class T>
void CkOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    #ifdef DISTRIBUTED_BASIS
 
    #else
    std::vector<pairIdx_t> finalIndList;
    std::vector<cdouble> factorList;
    this->Bf->genSymm(rowID, finalIndList, factorList);
    if (pm == LADDER::MINUS) {
        for (int i = 0; i < (int)finalIndList.size(); ++i){
            for (int r = 0; r < (int)posList.size(); ++r){
                cdouble factor = factorList.at(i) * expFactor.at(r);
                this->cp(spin, posList[r], factor, finalIndList[i], &rowMaps[0]);
            }
        }
    } else if (pm == LADDER::PLUS) {
        for (int i = 0; i < (int)finalIndList.size(); ++i){
            for (int r = 0; r < (int)posList.size(); ++r){
                cdouble factor = factorList.at(i) * expFactor.at(r);
                this->cm(spin, posList[r], factor, finalIndList[i], &rowMaps[0]);
            }
        } 
    }
    #endif
}

template<class T>
RamanOp<T>::RamanOp(Geometry* latt, Basis* Bi, int spmNum, int dmNum):OperatorBase<T>(latt, Bi, Bi, spmNum, dmNum) {
    pIn = VecD{1.0,0.0,0.0};
    pOut = VecD{0.0,1.0,0.0};
    NoRW = true;
}

template <class T>
void RamanOp<T>::setplz(VecD pIn_, VecD pOut_) {
    pIn = pIn_;
    pOut = pOut_;
    RamanWeight.resize(this->superExchangeJ.size());
    double norm = std::sqrt((this->latt->RdotR(pIn,pIn)) * (this->latt->RdotR(pOut,pOut)));
    int linkid = 0;
    for(auto& link : this->superExchangeJ){
        RamanWeight[linkid].resize(link.getLinkVecNum());
        for(int vecid=0; vecid < link.getLinkVecNum(); ++vecid){
            VecD r = link.getvec(vecid);
            RamanWeight[linkid][vecid] = (this->latt->RdotR(pIn,r))*(this->latt->RdotR(pOut,r))/(this->latt->RdotR(r,r))/norm;
        }
        ++linkid;
    }
    NoRW = false;
}
template <class T>
void RamanOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if (NoRW) setplz(pIn, pOut);
    std::vector<idx_t> finalIndList;
    std::vector<cdouble> factorList;
    this->pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < (int)finalIndList.size(); ++i) {
        int linkid = 0;
        for (const auto& link : this->superExchangeJ) {
            int matID = link.getmatid();
            int matIDp = matID; 
            if (link.isOrdered()) ++matIDp;
            T factor = factorList[i] * link.getVal();
            int bondid = 0;
            for (const auto& bond : link.bond()) {
                double rw = RamanWeight.at(linkid).at(link.getvecid(bondid));
                T factor_rw = factor * rw;
                int siteI = bond.at(0);
                int siteJ = bond.at(1);
                // szi * szj - 1/4*ni*nj 
                // szsznn(siteI, siteJ, factor_rw, finalIndList[i], &rowMaps[matID]);
                this->szsz(siteI, siteJ, factor_rw, finalIndList[i], &rowMaps[matID]);
                // 1/2 * sm.siteID * sp.siteIDp
                this->spsm(siteI, siteJ, factor_rw/2.0, finalIndList[i], &rowMaps[matID]);
                // 1/2 * sp.siteID * sm.siteIDp
                this->smsp(siteI, siteJ, factor_rw/2.0, finalIndList[i], &rowMaps[matID]);
                bondid++;
            }
            linkid++;
        }
    }
}

/*
    ********************
    * Correlator Class *
    ********************
*/
template <class T>
SSOp<T>::SSOp(Geometry* latt, Basis* Bi, int spmNum, int dmNum, int spindim):OperatorBase<T>(latt, Bi, Bi, spmNum, dmNum),\
    r(-1), siteJList(latt->getSiteNum()) {
    VecD coordi(3), coordr(3), coordf(3);
    for (int rIndex = 0; rIndex < latt->getSiteNum(); ++rIndex) {
        siteJList.at(rIndex).resize(latt->getSiteNum());
        latt->getSiteR(rIndex, coordr.data());
        for (int i = 0; i < latt->getOrbNum(); ++i){
            latt->getOrbR(i, coordi.data());
            vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), 3);
            int siteJ;
            if (latt->coordToOrbid(latt->getOrb(i), coordf.data(), siteJ)) {
                siteJList.at(rIndex).at(i) = siteJ;
            } else {
                std::cout<<"translation position not found for orbid = "<<i<<", transVecid = "<<rIndex<<"\n";
                exit(1);
            }
        }
    }
}

template <class T>
void SSOp<T>::setr(int r_) {
    assert(r_ < this->latt->getSiteNum());
    r=r_;
}

template <class T>
void SSOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    //binary rep
    if (this->Bf->getSiteDim()==2) {
        int kIndex = this->Bf->getkIndex();
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->Bf->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < (int)finalIndList.size(); ++i) {
            T factor = factorList.at(i);
            for (int siteI = 0; siteI < this->latt->getOrbNum(); ++siteI) {
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    this->szsz(siteI, siteJ, factor, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                } else {
                    for (int rIndex = 0; rIndex < this->latt->getSiteNum(); ++rIndex){
                        int siteJ = siteJList[rIndex][siteI];
                        // sz.siteI * sz.siteJ
                        this->szsz(siteI, siteJ, factor, finalIndList[i], &rowMaps[0]);
                        // 1/2 * sm.siteI * sp.siteJ
                        this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                        // 1/2 * sp.siteI * sm.siteJ
                        this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                    }
                }
            }
        }
    } else {
        int kIndex = this->Bf->getkIndex();
        VecI initVec(this->latt->getOrbNum());
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->Bf->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < (int)finalIndList.size(); ++i) {
            this->Bf->repToVec(finalIndList[i], initVec);
            T factor = factorList.at(i);
            for (int siteI = 0; siteI < this->latt->getOrbNum(); ++siteI) {
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    this->szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                } else {
                    for (int rIndex = 0; rIndex < this->getSiteNum(); ++rIndex) {
                        int siteJ = siteJList[rIndex][siteI];
                        // sz.siteI * sz.siteJ
                        this->szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMaps[0]);
                        // 1/2 * sm.siteI * sp.siteJ
                        this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                        // 1/2 * sp.siteI * sm.siteJ
                        this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                    }
                }
            }
        }
    }
}

template <class T>
void SSOp<T>::project(double s, T* vec){
    assert(r == -1);
    double stot = s * ( s + 1 );
    std::vector<T> tmp(this->getnloc());
    double smin = double(this->latt->getOrbNum() % 2) / 2.0;
    double smax = double(this->latt->getOrbNum()) / 2.0 + 0.1;
    for(double si = smin; si < smax; si += 1.0) {
        if (std::abs(si-s) > 0.1) {
            MxV(vec, tmp.data());
            double stoti = si * (si+1.0);
            #pragma omp parallel for
            for(idx_t i =0; i < this->getnloc(); ++i) {
                vec[i] = (tmp[i] - stoti * vec[i]) / (stot - stoti);
            }
        }
    }
}

template <class T>
SzkOp<T>::SzkOp(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum, int dmNum, int spindim):OperatorBase<T>(latt, Bi, Bf, spmNum, dmNum), expFactor(latt->getSiteNum()) {
    assert(Bi->getPGIndex()==-1 and Bf->getPGIndex()==-1);
    Ki = Bi->getkIndex();
    Kf = Bf->getkIndex();
    // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
    for (int i = 0; i < latt->getSiteNum(); ++i) {
        expFactor[i] = latt->expKR(Kf,i)/latt->expKR(Ki,i);
    }
}

template <class T>
void SzkOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if (this->Bi->getSiteDim() == 2) {
        #ifdef DISTRIBUTED_BASIS
            idx_t repI = this->Bf->getRepI(rowID);
            T dval = 0.0;
            for (int siteID = 0; siteID < this->latt->getOrbNum(); ++siteID) {
                dval += this->getSz(siteID, repI) * expFactor[siteID];
            }
            dval *= this->Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            idx_t colID;
            auto repI = this->Bf->getRepI(rowID);
            if (this->Bi->search(repI, colID)) {
                T dval = 0.0;
                for (int siteID = 0; siteID < this->latt->getOrbNum(); ++siteID) {
                    dval += this->getSz(siteID,repI) * expFactor[siteID];
                }
                dval *= this->Bf->getNorm(rowID)/this->Bi->getNorm(colID);
                rowMaps[0][colID] = dval;
            }
        #endif
    } else {
        #ifdef DISTRIBUTED_BASIS
            idx_t repI = this->Bf->getRepI(rowID);
            VecI initVec(this->latt->getOrbNum());
            T dval = 0.0;
            this->Bf->repToVec(repI, initVec);
            for (int siteID = 0; siteID < this->latt->getOrbNum(); ++siteID) {
                dval += this->getSz(siteID,initVec) * expFactor[siteID];
            }
            dval *= this->Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            idx_t colID;
            auto repI = this->Bf->getRepI(rowID);
            if (this->Bi->search(repI, colID)) {
                VecI initVec(this->latt->getOrbNum());
                T dval = 0.0;
                this->Bi->repToVec(repI, initVec);
                for (int siteID = 0; siteID < this->latt->getOrbNum(); ++siteID) {
                    dval += this->getSz(siteID,initVec) * expFactor[siteID];
                }
                dval *= this->Bf->getNorm(rowID)/this->Bi->getNorm(colID);
                rowMaps[0][colID] = dval;
            }
        #endif
    }
}

#endif // Operators_hpp