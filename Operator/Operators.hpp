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

/* Operators_hpp */

/***************
 * HAMILTONIAN *
 ***************/

template <LATTICE_MODEL MODEL, typename T> 
class Hamiltonian : public OperatorBase<T> {
public:
    Hamiltonian( ) { }
    Hamiltonian(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_ = 1, int dmNum_ = 0);
    ~Hamiltonian( ) { }

    // Add onsite energy V
    void pushV(std::vector<ORBITAL> orbList, double val);

    // Add onsite Coulomb interaction U
    void pushU(std::vector<ORBITAL> orbList, double val);

    void printV(std::ostream& os) const;

    void printU(std::ostream& os) const;

    // Calculate the sum of onsite energy and Coulomb interaction
    double diagVal(const VecI& occ, const VecI& docc) const;
    
    void setPeierls(Pulse* pulse = nullptr);

    void printPeierls(std::ostream& os = std::cout);

    bool next( );

    void row(idx_t rowidx, std::vector<MAP<T>>& rowMaps);

private:
    VecD V, U;
    Pulse* pulse{nullptr};
    VecD PeierlsOverlap;
};


class Current: public OperatorBase<cdouble> {
public:
    Current(Geometry *latt, Basis *Ba);
    ~Current(){};

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
class RamanOp: public SpinOperator<T>, public SparseMatrix<T>{
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
        SpinOperator<T>(pt_Ba),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_), pt_lattice(pt_lat),linkCount(0),spmCount(0){
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
        for (int i = 0; i < Links_.size(); ++i) pushLink(Links_[i], spmCount);
        spmCount++;
        assert(spmCount<=SparseMatrix<T>::spmNum);
        return *this;
    }

    void setplz(VecD pIn_, VecD pOut_);
    void setVal(int matID, double val){(SparseMatrix<T>::parameters).at(matID) = val;}
    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);
};

template <class T>
class SzkOp: public SpinOperator<T>, public SparseMatrix<cdouble>{
private:
    // initial k index Ki and final k index Kf
    int Ki, Kf;
    Basis *pt_Bi,*pt_Bf;
    Geometry *pt_lattice;
    std::vector<cdouble> expFactor;
public:
    SzkOp(Geometry *pt_lat, Basis *pt_Bi_, Basis *pt_Bf_, int spmNum_=1, int spindim=2):pt_Bi(pt_Bi_),pt_Bf(pt_Bf_), pt_lattice(pt_lat),expFactor(pt_lattice->getSiteNum()),\
        SpinOperator<T>(pt_Bi_),SparseMatrix<cdouble>(pt_Bi_,pt_Bf_,pt_Bf_->getSubDim(),spmNum_){
            assert(pt_Bi->getPGIndex()==-1 and pt_Bf->getPGIndex()==-1);
            Ki = pt_Bi->getkIndex();
            Kf = pt_Bf->getkIndex();
            // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
            for (int i = 0; i < pt_lattice->getSiteNum(); ++i) {
                expFactor[i] = pt_lattice->expKR(Kf,i)/pt_lattice->expKR(Ki,i);
            }
        };
    ~SzkOp(){};
    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);
    // void genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q);
    void genMat();
};

// correlator:Si*Si+r
// can also be used to construct total spin operator: S^2
template <class T>
class SSOp: public SpinOperator<T>, public SparseMatrix<T>{
private:
    Geometry *pt_lattice;
    int r; // if r = -1, sum.r sum.i Si*Si+r. total S^2
    std::vector<VecI> siteJList;
public:
    SSOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_=1, int spindim=2);
    ~SSOp(){}
    void setr(int r_){assert(r_<pt_lattice->getSiteNum());r=r_;}
    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);
    // S(i)*S(i+r)
    void genPairMat(int rIndex);
    void project(double s, T* vec);
};

/*
    ****************************
    * Operator Implementations *
    ****************************
*/

template<LATTICE_MODEL MODEL, typename T>
Hamiltonian<MODEL, T>::Hamiltonian(Geometry* latt, Basis* Bi, Basis* Bf, int spmNum_, int dmNum_):\
OperatorBase<T>(latt, Bi, Bf, spmNum_, dmNum_), V(latt->getUnitOrbNum(),0.0), U(latt->getUnitOrbNum(),0.0) {

}

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::pushV(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbID(orb);
        for (auto& id : ids) {
            V.at(id) = val;
        }
    }
}

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::pushU(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbID(orb);
        for (auto& id : ids) {
            U.at(id) = val;
        }
    }
}

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::setPeierls(Pulse* pulse) {
    if constexpr (MODEL == LATTICE_MODEL::HUBBARD or MODEL == LATTICE_MODEL::tJ) {
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
            for (int idx = 0; idx < overlap.size(); ++idx) {
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
        std::cout<<"Peierls substitution is not defined for "<<MODEL<<"\n";
        return;
    }
}

template<LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::printPeierls(std::ostream& os) {
    if constexpr (MODEL == LATTICE_MODEL::HUBBARD or MODEL == LATTICE_MODEL::tJ) {
        for (size_t id = 0; id < PeierlsOverlap.size(); ++id) {
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


template<LATTICE_MODEL MODEL, typename T>
bool Hamiltonian<MODEL, T>::next( ) {
    if (pulse) {
        auto A = pulse->getAa();
        for (int i = 1; i < PeierlsOverlap.size(); ++i) {
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

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::printV(std::ostream& os) const {
    os<<"V: ";
    for (auto val : V) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::printU(std::ostream& os) const {
    os<<"U: ";
    for (auto val : U) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <LATTICE_MODEL MODEL, typename T>
double Hamiltonian<MODEL, T>::diagVal(const VecI& occ, const VecI& docc) const {
    double val{0.0};
    for (int i = 0; i < this->latt->getUnitOrbNum(); ++i) {
        val += occ.at(i) * V.at(i) + docc.at(i) * U.at(i);
    }
    return val;
}

template <LATTICE_MODEL MODEL, typename T>
void Hamiltonian<MODEL, T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if constexpr (MODEL == LATTICE_MODEL::HUBBARD) {
        // diagonal part. occupancy and double-occ
        VecI occ, docc;
        idx_t repI = this->Bf->getRepI(rowID);
        pairIdx_t pairRepI = this->Bf->getPairRepI(repI);
        this->latt->orbOCC(pairRepI, occ, docc);
        T val = diagVal(occ,docc);
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
    } else if constexpr(MODEL == LATTICE_MODEL::tJ) {
        // off diagonal part
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->Bf->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); ++i){
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
    } else if constexpr(MODEL == LATTICE_MODEL::HEISENBERG) {
        if(this->Bf->getSiteDim()==2){
            int kIndex = this->Bf->getkIndex();
            std::vector<idx_t> finalIndList;
            std::vector<cdouble> factorList;
            this->Bf->genSymm(rowID, finalIndList, factorList);
            for (int i = 0; i < finalIndList.size(); ++i){
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
            int kIndex = this->Bf->getkIndex();
            VecI initVec(this->latt->getOrbNum());
            std::vector<idx_t> finalIndList;
            std::vector<cdouble> factorList;
            this->Bf->genSymm(rowID, finalIndList, factorList);
            for (int i = 0; i < finalIndList.size(); ++i){
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
    assert(posList.size() == this->latt->getSiteNum());
}

template <class T>
void CkOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    #ifdef DISTRIBUTED_BASIS
 
    #else
    std::vector<pairIdx_t> finalIndList;
    std::vector<cdouble> factorList;
    this->Bf->genSymm(rowID, finalIndList, factorList);
    if (pm == LADDER::MINUS) {
        for (int i = 0; i < finalIndList.size(); ++i){
            for (int r = 0; r < posList.size(); ++r){
                cdouble factor = factorList.at(i) * expFactor.at(r);
                this->cp(spin, posList[r], factor, finalIndList[i], &rowMaps[0]);
            }
        }
    } else if (pm == LADDER::PLUS) {
        for (int i = 0; i < finalIndList.size(); ++i){
            for (int r = 0; r < posList.size(); ++r){
                cdouble factor = factorList.at(i) * expFactor.at(r);
                this->cm(spin, posList[r], factor, finalIndList[i], &rowMaps[0]);
            }
        } 
    }
    #endif
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
void RamanOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    if(NoRW) setplz(pIn, pOut);
    std::vector<idx_t> finalIndList;
    std::vector<cdouble> factorList;
    this->pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); ++i){
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
SSOp<T>::SSOp(Geometry *pt_lat, Basis *pt_Ba, int spmNum_, int spindim):r(-1),pt_lattice(pt_lat), siteJList(pt_lat->getSiteNum()),\
    SpinOperator<T>(pt_Ba),SparseMatrix<T>(pt_Ba,pt_Ba,pt_Ba->getSubDim(),spmNum_){
    VecD coordi(3), coordr(3), coordf(3);
    for (int rIndex = 0; rIndex < pt_lat->getSiteNum();rIndex++){
        siteJList.at(rIndex).resize(pt_lat->getSiteNum());
        pt_lattice->getSiteR(rIndex, coordr.data());
        for (int i = 0; i < pt_lattice->getOrbNum(); ++i){
            pt_lattice->getOrbR(i,coordi.data());
            vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), 3);
            int siteJ;
            if (pt_lattice->coordToOrbid(pt_lattice->getOrb(i), coordf.data(), siteJ)){
                siteJList.at(rIndex).at(i) = siteJ;
            }else{
                std::cout<<"translation position not found for orbid = "<<i<<", transVecid = "<<rIndex<<std::endl;
                exit(1);
            }
        }
    }
}

template <class T>
void SSOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    //binary rep
    if(this->pt_Basis->getSiteDim()==2){
        int kIndex = this->pt_Basis->getkIndex();
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); ++i){
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); ++siteI){
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    this->szsz(siteI, siteJ, factor, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], &rowMaps[0]);
                }
                else{
                    for (int rIndex = 0; rIndex < pt_lattice->getSiteNum(); rIndex++){
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
    }
    else{
        int kIndex = this->pt_Basis->getkIndex();
        VecI initVec(pt_lattice->getOrbNum());
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); ++i){
            this->pt_Basis->repToVec(finalIndList[i], initVec);
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); ++siteI){
                if (r >= 0){
                    int siteJ = siteJList[r][siteI];
                    // sz.siteI * sz.siteJ
                    this->szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sm.siteI * sp.siteJ
                    this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                    // 1/2 * sp.siteI * sm.siteJ
                    this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMaps[0]);
                }
                else{
                    for (int rIndex = 0; rIndex < pt_lattice->getSiteNum(); rIndex++){
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
    assert(r==-1);
    double stot = s*(s+1);
    std::vector<T> tmp(BaseMatrix<T>::getnloc());
    double smin = double(pt_lattice->getOrbNum()%2)/2.0;
    double smax = double(pt_lattice->getOrbNum())/2.0+0.1;
    for(double si = smin; si<smax; si += 1.0){
        if (std::abs(si-s)>0.1){
            MxV(vec, tmp.data());
            double stoti = si * (si+1.0);
            #pragma omp parallel for
            for(idx_t i =0; i < BaseMatrix<T>::getnloc();++i) vec[i] = (tmp[i]-stoti*vec[i])/(stot-stoti);
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
    int kIndex = this->pt_Basis->getkIndex();
    this->clear();
    this->reserve(pt_lattice->getOrbNum()/2+1);
    MAP<T> rowMap;
    this->pushRow(&rowMap);
    VecI initVec(pt_lattice->getOrbNum());
    for (idx_t rowID = this->startRow; rowID < this->endRow; rowID++){
        rowMap.clear();
        std::vector<idx_t> finalIndList;
        std::vector<cdouble> factorList;
        this->pt_Basis->genSymm(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); ++i){
            this->pt_Basis->repToVec(finalIndList[i], initVec);
            cdouble factor = factorList.at(i);
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); ++siteI){
                int siteJ = siteJList[rIndex][siteI];
                // sz.siteI * sz.siteJ
                this->szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMap);
                // 1/2 * sm.siteI * sp.siteJ
                this->spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
                // 1/2 * sp.siteI * sm.siteJ
                this->smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
            }
        }
        this->pushRow(&rowMap);
    }
    #endif
}

template <class T>
void SzkOp<T>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    if(pt_Bi->getSiteDim()==2){
        #ifdef DISTRIBUTED_BASIS
            idx_t repI = pt_Bf->getRepI(rowID);
            cdouble dval = 0.0;
            for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                dval += getSz(siteID,repI) * expFactor[siteID];
            }
            dval *= pt_Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            idx_t colID;
            if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                // std::cout<<"colID:"<<colID<<"\n";
                idx_t repI = pt_Bi->getRepI(colID);
                cdouble dval = 0.0;
                for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                    dval += this->getSz(siteID,repI) * expFactor[siteID];
                }
                dval *= pt_Bf->getNorm(rowID)/pt_Bi->getNorm(colID);
                // std::cout<<"rowMaps size:"<<rowMaps.size()<<"\n";
                rowMaps[0][colID] = dval;
            }
        #endif
    }
    else{
        #ifdef DISTRIBUTED_BASIS
            idx_t repI = pt_Bf->getRepI(rowID);
            VecI initVec(pt_lattice->getOrbNum());
            cdouble dval = 0.0;
            pt_Bf->repToVec(repI, initVec);
            for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                dval += getSz(siteID,initVec) * expFactor[siteID];
            }
            dval *= pt_Bf->getNorm(rowID);
            rowMaps[0][repI] = dval;
        #else
            idx_t colID;
            if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                VecI initVec(pt_lattice->getOrbNum());
                cdouble dval = 0.0;
                pt_Bi->repToVec(pt_Bi->getRepI(colID), initVec);
                for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                    dval += this->getSz(siteID,initVec) * expFactor[siteID];
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
    MAP<T> rowMap;
    pushRow(&rowMap);
    cdouble dval;
    VecI initVec(pt_lattice->getOrbNum());
    switch(PARTITION){
        case ROW_PARTITION:{
            idx_t colID;
            for (idx_t rowID = startRow; rowID < endRow; rowID++){
                rowMap.clear();
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bi->repToVec(pt_Bi->getRepI(colID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += this->getSz(siteID,initVec) * expFactor[siteID];
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
            idx_t colID;
            for (idx_t rowID = startRow; rowID < endRow; rowID++){
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bf->repToVec(pt_Bf->getRepI(rowID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += this->getSz(siteID,initVec) * expFactor[siteID];
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