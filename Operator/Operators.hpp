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
#include "Pulse/pulse.hpp"
#include "Solver/PARPACKSolver.hpp"

/***************
 * HAMILTONIAN *
 ***************/

template <typename T, IsBasisType B>
class HamiltonianBase : public OperatorBase<T, B> {
public:
    HamiltonianBase( ) = default;
    HamiltonianBase(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool commuteWithTrans = true, bool commuteWithPG = true, int spmNum_ = 1, int dmNum_ = 0);
    ~HamiltonianBase( ) = default;

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

    virtual void row(idx_t rowidx, std::vector<MAP<T>>& rowMaps){};

    void diag(int nev = 1);
    T getEval(int n = 0);
    T* getEvec(int n = 0);

protected:
    PARPACKSolver<T> solver;
    VecD V, U;
};

template <typename T, IsBasisType B>
class Hamiltonian: public HamiltonianBase<T, B> {
public:
    Hamiltonian(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool commuteWithTrans = true, bool commuteWithPG = true, int spmNum_ = 1, int dmNum_ = 1);
    ~Hamiltonian( ) = default;

    void setPeierls(Pulse* pulse = nullptr);

    void printPeierls(std::ostream& os = std::cout);

    bool next( );

    void row(idx_t rowidx, std::vector<MAP<T>>& rowMaps);

private:
    Pulse* pulse{nullptr};
    VecD PeierlsOverlap;
};


template<ContainCharge B>
class Current: public OperatorBase<cdouble, B> {
public:
    Current(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool commuteWithTrans = true, bool commuteWithPG = false);

    void setDirection(const std::string &plz);

    void print(std::ostream &os = std::cout) const;

    void row(idx_t rowID, std::vector<MAP<cdouble>>& rowMaps);

private:
    void setDirection(Vec3d d);

    Vec3d direction;
    std::vector<Link<cdouble>> myHoppingT;
    std::string plz;
};

template<ContainCharge B>
class Nocc: public OperatorBase<double, B> {
public:
    Nocc(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool commuteWithTrans = true, bool commuteWithPG = true);
    ~Nocc( ) { }
    
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
    void save(const std::string &dir);

private:
    std::vector<std::vector<double>> records;
};

template <class T, ContainCharge B>
class CkOp: public OperatorBase<T, B> {
public:
    CkOp(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool commuteWithTrans = false, bool commuteWithPG = false);

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

template<class T, IsPureSpin B>
class RamanOp: public OperatorBase<T, B> {
public:
    RamanOp(Geometry *pt_lat, Basis<B> *Bi, Basis<B> *Bf, bool commuteWithTrans = true, bool commuteWithPG = false, int spmNum_=1, int dmNum = 0);

    void setplz(Vec3d pIn_, Vec3d pOut_);

    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);

private:
    // polarization of in/out photons
    Vec3d pIn, pOut;
    std::vector<VecD> RamanWeight;
    bool NoRW;
};

// correlator:Si*Si+r
// can also be used to construct total spin operator: S^2
template <class T, ContainSpin B>
class SSOp: public OperatorBase<T, B> {
public:
    SSOp(Geometry *latt, Basis<B>* Bi, Basis<B> *Bf, bool commuteWithTrans = true, bool commuteWithPG = false, int spmNum = 1, int dmNum = 0);

    void setr(int r_);

    void row(idx_t rowID, std::vector<MAP<T>>& rowMaps);
 
    void project(double s, T* vec);

private:
    // if r = -1, sum.r sum.i Si*Si+r. total S^2
    int r; 
    std::vector<VecI> siteJList;
};

template <class T, ContainSpin B>
class SzkOp: public OperatorBase<T, B> {
public:
    SzkOp(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool commuteWithTrans = false, bool commuteWithPG = false, int spmNum = 1, int dmNum = 0);

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

template<typename T, IsBasisType B>
HamiltonianBase<T, B>::HamiltonianBase(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool trans, bool pg, int spmNum_, int dmNum_):\
OperatorBase<T, B>(latt, Bi, Bf, trans, pg, spmNum_, dmNum_), V(latt->getOrbNum(),0.0), U(latt->getOrbNum(),0.0) {

}

template<typename T, IsBasisType B>
void HamiltonianBase<T, B>::diag(int nev) {
    solver = PARPACKSolver<T>(this, nev);
    solver.diag();
}

template<typename T, IsBasisType B>
T HamiltonianBase<T, B>::getEval(int n) {
    return solver.getEigval(n);
}

template<typename T, IsBasisType B>
T* HamiltonianBase<T, B>::getEvec(int n) {
    return solver.getEigvec(n);
}

template <typename T, IsBasisType B>
void HamiltonianBase<T, B>::pushV(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbPos(orb);
        for (auto& id : ids) {
            V.at(id) = val;
        }
    }
}

template <typename T, IsBasisType B>
void HamiltonianBase<T, B>::pushU(std::vector<ORBITAL> orbList, double val) {
    for (const auto& orb : orbList) {
        auto ids = this->latt->getOrbPos(orb);
        for (auto& id : ids) {
            U.at(id) = val;
        }
    }
}

template <typename T, IsBasisType B>
void HamiltonianBase<T, B>::printV(std::ostream& os) const {
    os<<"V: ";
    for (auto val : V) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <typename T, IsBasisType B>
void HamiltonianBase<T, B>::printU(std::ostream& os) const {
    os<<"U: ";
    for (auto val : U) {
        os<<val<<", ";
    }
    os<<"\n";
}

template <typename T, IsBasisType B>
double HamiltonianBase<T, B>::diagVal(const VecI& occ, const VecI& docc) const {
    double val{0.0};
    for (int i = 0; i < this->latt->getUnitOrbNum(); ++i) {
        val += occ.at(i) * V.at(i) + docc.at(i) * U.at(i);
    }
    return val;
}

template <typename T, IsBasisType B>
Hamiltonian<T, B>::Hamiltonian(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf,  bool trans, bool pg, int spmNum, int dmNum) : HamiltonianBase<T,B>(latt, Bi, Bf, trans, pg, spmNum, dmNum) {
    
}

template <typename T, IsBasisType B>
void Hamiltonian<T, B>::setPeierls(Pulse* pulse) {
    if constexpr (ContainCharge<B>) {
        this->pulse = pulse;
        if (pulse) {
            auto epol = pulse->getPol();
            normalize(epol);
            VecD overlap;
            for(const auto& link : this->hoppingT) {
                assert_msg(link.getLinkVecNum()==1, "In setPeierls, each hopping link should only have one dr!");
                auto dr = link.getvec(0);
                dr = this->latt->RtoRxy(dr);
                auto val = dot(epol,dr);
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
        std::cout << "Peierls substitution is not defined for charge!\n";
        return;
    }
}

template<typename T, IsBasisType B>
void Hamiltonian<T, B>::printPeierls(std::ostream& os) {
    if constexpr (ContainCharge<B>) {
        for (int id = 0; id < (int)PeierlsOverlap.size(); ++id) {
            os << "A * dr = " << PeierlsOverlap[id] << std::endl;
            auto matid = (id == 0) ? id : 2 * id -1;
            for (const auto& link : this->hoppingT) {
                if (link.getmatid() == matid) {
                    link.print();
                }
            }
            os << std::endl;
        }
    }
}

template<typename T, IsBasisType B>
bool Hamiltonian<T, B>::next( ) {
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

template <typename T, IsBasisType B>
void Hamiltonian<T, B>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if constexpr (ContainCharge<B>) {
        auto state = this->Bf->get(rowID);
        auto nf = this->Bf->norm(rowID);
        // diagonal part. occupancy and double-occ
        T val = 0.0;
        for (int i = 0; i < B::getNSite(); ++i) {
            if (!this->V.empty()) {
                if (bitTest(state.upState(), i)) {
                    val += this->V.at(i);
                }
                if (bitTest(state.upState(), i)) {
                    val += this->V.at(i);
                }
            }
            if (!this->U.empty()) {
                if (bitTest(state.upState(), i) && bitTest(state.dnState(), i)) {
                    val += this->U.at(i);
                }
            }
        }
        for (const auto& link : this->interBandU){
            for (const auto& bond : link.bond()){
                auto siteI = bond.at(0);
                auto siteJ = bond.at(1);
                assert(siteI!=siteJ);
                val += link.getVal()
                        * double(bitTest(state.upState(), siteI) + bitTest(state.dnState(), siteI))
                        * double(bitTest(state.upState(), siteJ) + bitTest(state.dnState(), siteJ));
            }
        }
        SparseMatrix<T>::putDiag(val,rowID);
        // off diagonal part
//        for (const auto& trOp : this->trHoppingT) {
//            auto statef = state;
//            statef.transform(trOp.g);
//            for (const auto& bond : trOp.Op.bonds) {
//                auto val = bond.val/nf;
//                this->pushElement( val * (CPlus<SPIN::UP>(bond[0]) * (CMinus<SPIN::UP>(bond[1]) * BVopt<T, B>(statef))), mapPtr );
//            }
//        }
        for (const auto& link:this->hoppingT) {
            int matID = link.getmatid();
            int matIDp = matID;
            if (link.isOrdered()) matIDp++;
            T val = link.getVal() / nf;
            for (const auto& bond : link.bond()) {
                int siteI = bond.at(0);
                int siteJ = bond.at(1);
                cdouble phase = this->latt->twistPhase(siteI,siteJ);
                cdouble phase_c = std::conj(phase);
                // cp.siteI * cm.siteJ
                this->pushElement( phase * val * (CPlus<SPIN::UP>(siteI) * (CMinus<SPIN::UP>(siteJ) * BVopt<T, B>(state))),  &rowMaps[matID]);
                this->pushElement( phase_c * val * (CPlus<SPIN::UP>(siteJ) * (CMinus<SPIN::UP>(siteI) * BVopt<T, B>(state))),  &rowMaps[matIDp]);
                this->pushElement( phase * val * (CPlus<SPIN::DOWN>(siteI) * (CMinus<SPIN::DOWN>(siteJ) * BVopt<T, B>(state))),  &rowMaps[matID]);
                this->pushElement( phase_c * val * (CPlus<SPIN::DOWN>(siteJ) * (CMinus<SPIN::DOWN>(siteI) * BVopt<T, B>(state))),  &rowMaps[matIDp]);
            }
        }
//
//        for (const auto& link:this->exchangeJ) {
//            int matID = link.getmatid();
//            cdouble factor = factorList.at(i) * link.getVal();
//            for (const auto& bond : link.bond()){
//                int siteI = bond.at(0);
//                int siteJ = bond.at(1);
//                this->exchange(siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
//            }
//        }
//
//        for (const auto& link:this->pairHoppingJ) {
//            int matID = link.getmatid();
//            cdouble factor = factorList.at(i) * link.getVal();
//            for (const auto& bond : link.bond()){
//                int siteI = bond.at(0);
//                int siteJ = bond.at(1);
//                this->pairHopping(siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
//            }
//        }
    }
    if constexpr(ContainSpin<B>) {
        int matID = 0;
        MAP<T>* mapPtr = &rowMaps[matID];
        auto repI0 = this->Bf->get(rowID);
        auto nf = this->Bf->norm(rowID);
        std::cout << "init state: " << repI0 << ", norm:" << nf;
        for (const auto& trOp : this->trSuperExchangeJ) {
            auto repI = repI0;
            repI.transform(trOp.g);
            for (const auto& bond : trOp.Op.bonds) {
                std::cout << "bond: " << bond.val << " * " << bond.sites[0] << "--" << bond.sites[1] << std::endl;
                auto val = bond.val/nf;
                this->pushElement( val * (Sz(bond[0]) * (Sz(bond[1]) * BVopt<T, B>(repI))), mapPtr );
                this->pushElement( val / 2.0 * (SPlus(bond[0]) * (SMinus(bond[1]) * BVopt<T, B>(repI))), mapPtr );
                this->pushElement( val / 2.0 * (SMinus(bond[0]) * (SPlus(bond[1]) * BVopt<T, B>(repI))), mapPtr );
            }
        }
        //TODO:Chiral Terms
//        for (const auto& trOp : this->trChiralTermK) {
//            auto repI = repI0;
//            repI.transform(trOp.g);
//            for (const auto& bond : trOp.Op.bonds) {
//                this->chiral(bond[0], bond[1], bond[2], bond.val/nf, repI, mapPtr);
//            }
//        }
    }
    if constexpr (ContainPhonon<B>) {

    }
}

template<ContainCharge B>
Nocc<B>::Nocc(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool trans, bool pg) : OperatorBase<double, B>(latt, Bi, Bf, trans, pg, 0, latt->getUnitOrbNum()) {
    for(int i = 0; i < latt->getUnitOrbNum(); ++i)  {
        this->diagValList.resize(BaseMatrix<double>::nloc);
    }
    records.resize(latt->getUnitOrbNum());
}

template<ContainCharge B>
void Nocc<B>::row(idx_t rowID) {
    // diagonal part. occupancy
    VecI occ;
    auto state = this->Bf->get(rowID);
    this->latt->orbOCC({state.upState(), state.dnState()}, occ);
    idx_t loc_rowID = rowID - BaseMatrix<double>::startRow;
    for (int i = 0; i < this->latt->getUnitOrbNum(); ++i) {
        SparseMatrix<double>::diagValList[i][loc_rowID] = occ[i];
    }
}

template<ContainCharge B>
void Nocc<B>::construct( ) {
    #pragma omp parallel for
    for(idx_t rowID = BaseMatrix<double>::startRow; rowID < BaseMatrix<double>::endRow; rowID++) {
        row(rowID);
    }
}

template<ContainCharge B>
double Nocc<B>::count(ORBITAL orbital, dataType* vec) {
    VecI ids = this->latt->getOrbID(orbital);
    double sum_final = 0.0;
    for (auto id : ids) {
        double sum = 0.0, part_sum = 0.0;
        #pragma omp parallel for reduction(+:part_sum)
        for(idx_t i = 0; i < BaseMatrix<double>::nloc; ++i){
            part_sum += this->diagValList[id][i] * (std::real(vec[i]) * std::real(vec[i]) + std::imag(vec[i]) * std::imag(vec[i]));
        }
        MPI_Allreduce(&part_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sum_final += sum;
    }
    return sum_final;
}

template<ContainCharge B>
double Nocc<B>::count(int id, dataType* vec) {
    double sum = 0.0;
    double part_sum = 0.0;
    #pragma omp parallel for reduction(+:part_sum)
    for(idx_t i = 0; i < BaseMatrix<double>::nloc; ++i){
        part_sum += this->diagValList[id][i] * (std::real(vec[i]) * std::real(vec[i]) + std::imag(vec[i]) * std::imag(vec[i]));
    }
    MPI_Allreduce(&part_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

template<ContainCharge B>
void Nocc<B>::count(dataType* vec) {
    for (int id = 0; id < this->latt->getUnitOrbNum(); ++id) {
        records.at(id).push_back(count(id, vec));
    }
}

template<ContainCharge B>
void Nocc<B>::clear( ) {
    records.clear();
}

template<ContainCharge B>
void Nocc<B>::save(const std::string &dir) {
    for (int id = 0; id < this->latt->getUnitOrbNum(); ++id) {
        ::save(records.at(id).data(), (int)records.at(id).size(), dir + "/orb" + tostr(id));
    }
}


template <class T, ContainCharge B>
CkOp<T, B>::CkOp(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool trans, bool pg):\
OperatorBase<T, B>(latt, Bi, Bf, trans, pg), expFactor(latt->getSiteNum()) {
    assert(Bi->getPGIndex()==-1 and Bf->getPGIndex()==-1);
    Ki = Bi->getkIndex();
    Kf = Bf->getkIndex();
    // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
    int N = latt->getSiteNum();
    for (int i = 0; i < N; ++i) {
        expFactor[i] = latt->expKR(Ki, i) / latt->expKR(Kf, i) / std::sqrt(N);
    }
}

template <class T, ContainCharge B>
void CkOp<T, B>::set(LADDER pm, SPIN spin, Orbital orb) {
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

//TODO: Check Ck with symmetry
template <class T, ContainCharge B>
void CkOp<T, B>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
    #ifdef DISTRIBUTED_BASIS
 
    #else
    auto state = this->Bf->get(rowID);
    auto nf = this->Bf->norm(rowID);
    if (pm == LADDER::MINUS) {
        for (int r = 0; r < (int)posList.size(); ++r){
            T factor = expFactor.at(r) / nf;
            if (spin == SPIN::UP) {
                this->pushElement(factor * (CPlus<SPIN::UP>(posList[r]) * BVopt<T, B>(state)), &rowMaps[0]);
            } else {
                this->pushElement(factor * (CPlus<SPIN::DOWN>(posList[r]) * BVopt<T, B>(state)), &rowMaps[0]);
            }
        }
    } else if (pm == LADDER::PLUS) {
        for (int r = 0; r < (int)posList.size(); ++r){
            T factor = expFactor.at(r) / nf;
            if (spin == SPIN::UP) {
                this->pushElement(factor * (CMinus<SPIN::UP>(posList[r]) * BVopt<T, B>(state)), &rowMaps[0]);
            } else {
                this->pushElement(factor * (CMinus<SPIN::DOWN>(posList[r]) * BVopt<T, B>(state)), &rowMaps[0]);
            }
        }
    }
    #endif
}

template<ContainCharge B>
Current<B>::Current(Geometry *latt, Basis<B> *Bi, Basis<B> *Bf, bool trans, bool pg) : OperatorBase<cdouble, B>(latt, Bi, Bf, trans, pg) {

}

template<ContainCharge B>
void Current<B>::setDirection(const std::string& plz)
{
    this->plz = plz;
    if (plz == "x")
        direction = {1.0, 0.0, 0.0};
    else if (plz == "y")
        direction = {0.0, 1.0, 0.0};
    else if (plz == "z")
        direction = {0.0, 0.0, 1.0};
    else
        direction = {0.0, 0.0, 0.0};
    setDirection(direction);
}

template<ContainCharge B>
void Current<B>::setDirection(Vec3d d) {
    direction = d;
    myHoppingT.clear();
    for (auto link : this->hoppingT) {
        assert_msg(link.getLinkVecNum()==1, "In setDirection, each hopping link should only have one dr!");
        auto dr = link.getvec(0);
        dr = this->latt->RtoRxy(dr);
        auto overlap = dot(direction, dr);
        if (overlap != 0.0) {
            link.setVal(link.getVal() * overlap);
            myHoppingT.push_back(link);
        }
    }
}

template<ContainCharge B>
void Current<B>::print(std::ostream &os) const {
    for (const auto& link : this->myHoppingT) {
        link.print(false, os);
    }
}

template<ContainCharge B>
void Current<B>::row(idx_t rowID, std::vector<MAP<cdouble>>& rowMaps) {
    auto state = this->Bf->get(rowID);
    auto nf = this->Bf->norm(rowID);
    for (const auto& link : this->myHoppingT) {
        int matID = link.getmatid();
        cdouble val = CPLX_I * link.getVal() / nf;
        for (const auto& bond : link.bond()){
            int siteI = bond.at(0);
            int siteJ = bond.at(1);
            // cp.siteI * cm.siteJ
            this->pushElement( val * (CPlus<SPIN::UP>(siteI) * (CMinus<SPIN::UP>(siteJ) * BVopt<cdouble, B>(state))),  &rowMaps[matID]);
            this->pushElement( -val * (CPlus<SPIN::UP>(siteJ) * (CMinus<SPIN::UP>(siteI) * BVopt<cdouble, B>(state))),  &rowMaps[matID]);
            this->pushElement( val * (CPlus<SPIN::DOWN>(siteI) * (CMinus<SPIN::DOWN>(siteJ) * BVopt<cdouble, B>(state))),  &rowMaps[matID]);
            this->pushElement( -val * (CPlus<SPIN::DOWN>(siteJ) * (CMinus<SPIN::DOWN>(siteI) * BVopt<cdouble, B>(state))),  &rowMaps[matID]);
        }
    }
}

template<class T, IsPureSpin B>
RamanOp<T, B>::RamanOp(Geometry* latt, Basis<B> *Bi, Basis<B> *Bf, bool trans, bool pg, int spmNum, int dmNum) : OperatorBase<T, B>(latt, Bi, Bf, trans, pg, spmNum, dmNum) {
    pIn = Vec3d{1.0,0.0,0.0};
    pOut = Vec3d{0.0,1.0,0.0};
    NoRW = true;
}

template<class T, IsPureSpin B>
void RamanOp<T, B>::setplz(Vec3d pIn_, Vec3d pOut_) {
    pIn = pIn_;
    pOut = pOut_;
    RamanWeight.resize(this->superExchangeJ.size());
    double norm = std::sqrt((this->latt->RdotR(pIn,pIn)) * (this->latt->RdotR(pOut,pOut)));
    int linkid = 0;
    for(auto& link : this->superExchangeJ){
        RamanWeight[linkid].resize(link.getLinkVecNum());
        for(int vecid=0; vecid < link.getLinkVecNum(); ++vecid){
            auto r = link.getvec(vecid);
            RamanWeight[linkid][vecid] = (this->latt->RdotR(pIn,r))*(this->latt->RdotR(pOut,r))/norm;
        }
        ++linkid;
    }
    NoRW = false;
}
template<class T, IsPureSpin B>
void RamanOp<T, B>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    if (NoRW) setplz(pIn, pOut);
    std::vector<idx_t> finalIndList;
    std::vector<cdouble> factorList;
    this->Bf->genSymm(rowID, finalIndList, factorList);
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
template <class T, ContainSpin B>
SSOp<T, B>::SSOp(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool trans, bool pg, int spmNum, int dmNum) : OperatorBase<T, B>(latt, Bi, Bf, trans, pg, spmNum, dmNum),\
    r(-1), siteJList(latt->getSiteNum()) {
    Vec3d coordi, coordr, coordf;
    for (int rIndex = 0; rIndex < latt->getSiteNum(); ++rIndex) {
        siteJList.at(rIndex).resize(latt->getSiteNum());
        coordr = latt->getSiteR(rIndex);
        for (int i = 0; i < latt->getOrbNum(); ++i){
            coordi = latt->getOrbR(i);
            coordf = coordi + coordr;
            int siteJ;
            if (latt->coordToOrbid(latt->getOrb(i), coordf, siteJ)) {
                siteJList.at(rIndex).at(i) = siteJ;
            } else {
                std::cout<<"translation position not found for orbid = "<<i<<", transVecid = "<<rIndex<<"\n";
                exit(1);
            }
        }
    }
}

template <class T, ContainSpin B>
void SSOp<T, B>::setr(int r_) {
    assert(r_ < this->latt->getSiteNum());
    r=r_;
}

//FIXME
template <class T, ContainSpin B>
void SSOp<T, B>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps){
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
}

//FIXME
template <class T, ContainSpin B>
void SSOp<T, B>::project(double s, T* vec){
    assert(r == -1);
    double stot = s * ( s + 1 );
    std::vector<T> tmp(this->getnloc());
    double smin = double(this->latt->getOrbNum() % 2) / 2.0;
    double smax = double(this->latt->getOrbNum()) / 2.0 + 0.1;
    for(double si = smin; si < smax; si += 1.0) {
        if (std::abs(si-s) > 0.1) {
            this->MxV(vec, tmp.data());
            double stoti = si * (si+1.0);
            #pragma omp parallel for
            for(idx_t i =0; i < this->getnloc(); ++i) {
                vec[i] = (tmp[i] - stoti * vec[i]) / (stot - stoti);
            }
        }
    }
}

template <class T, ContainSpin B>
SzkOp<T, B>::SzkOp(Geometry* latt, Basis<B>* Bi, Basis<B>* Bf, bool trans, bool pg, int spmNum, int dmNum) : OperatorBase<T, B>(latt, Bi, Bf, trans, pg, spmNum, dmNum), expFactor(latt->getSiteNum()) {
    // assert(Bi->getPGIndex()==-1 and Bf->getPGIndex()==-1);
    Ki = Bi->getkIndex();
    Kf = Bf->getkIndex();
    // expFactor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
    this->Sz.clear();
    //!Fix: this is Szk^\dagger 
    for (int i = 0; i < latt->getSiteNum(); ++i) {
        expFactor[i] = latt->expKR(Ki, i) / latt->expKR(Kf, i) / cdouble(latt->getSiteNum());
        this->Sz.add(Bond<T,1>(expFactor[i], {i}));
    }
    Generator<T> Gi, Gf;
    std::vector<Transform<T>> allTr;
    this->getGiGf(Gi, Gf, allTr);
    assignTrInteractions<T,1>(Gi, Gf, allTr, {this->Sz}, this->trSz, 'n');
}

//TODO: Check Szk with Symm
template <class T, ContainSpin B>
void SzkOp<T, B>::row(idx_t rowID, std::vector<MAP<T>>& rowMaps) {
    #ifdef DISTRIBUTED_BASIS
        idx_t repI = this->Bf->getRepI(rowID);
        T dval = 0.0;
        for (int siteID = 0; siteID < this->latt->getOrbNum(); ++siteID) {
            dval += this->getSz(siteID, repI) * expFactor[siteID];
        }
        dval *= this->Bf->getNorm(rowID);
        rowMaps[0][repI] = dval;
    #else
        auto repI = this->Bf->get(rowID);
        auto nf = this->Bf->norm(rowID);
        for (const auto& gSz : this->trSz) {
            auto repIf = repI;
            repIf.transform(gSz.g);
            auto colID =this->Bi->search(repIf);
            if (colID) {
                T val = 0.0;
                for (const auto& bond : gSz.Op.bonds) {
                    val += bond.val * this->getSz(bond[0], repIf);
                }
                val /= (nf * this->Bi->norm(colID));
                MapPush(&rowMaps.at(0), *colID, val);
            }
        }
    #endif
}

#endif // Operators_hpp