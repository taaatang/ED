#pragma once

#include <vector>
#include <cmath>
#include <iostream>

#include "basis/basis.hpp"

// push data to an unordered map
template <typename T>
inline void MapPush(MAP<T>* map_pt, idx_t key, T val){
    auto it = map_pt->find(key);
    if constexpr (std::is_same<cdouble, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = std::conj(val);
        } else {
            it->second += std::conj(val);
        }
    } else if constexpr (std::is_same<double, T>::value) {
        if (it == map_pt->end()) {
            (*map_pt)[key] = val;
        } else {
            it->second += val;
        }
    } else {
        std::cout<<"MapPush only defined for double/cdouble!\n";
        exit(1);
    }
}

template<typename T>
class FermionOperator{
public:
    FermionOperator( ) { }
    FermionOperator(Basis* pt_Ba, bool commuteWithTrans = false, bool commuteWithPG = false);
    ~FermionOperator( ) { };
    // c^dagger/c with spin act on siteI of pairRepI 
    bool cp(SPIN spin, int siteI, pairIdx_t& pairRepI, int &sign);
    bool cm(SPIN spin, int siteI,  pairIdx_t& pairRepI, int &sign);
    // hopping deals with sign count difference for hubbard and tJ
    bool cpcm(SPIN spin, int siteI, int siteJ, pairIdx_t& pairRepI, int& sign);

    // used for single particle spectra 
    void cp(SPIN spin, int siteI, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);
    void cm(SPIN spin, int siteI, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);
    // hopping term
    void cpcm(SPIN spin, int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);
   
    void diag(idx_t rowID, T factor, MAP<T>* rowMap);

    /**
     * @brief Coulomb interactions
     * https://doi.org/10.1016/S0370-1573(00)00121-6
     * i is the orbital index on the same site
     * Intraband Coulomb interaction: U * n_i_up * n_i_dn
     * Interband Coulomb interaction: U' * (n_i_up + n_i_dn) * (n_j_up + n_j_dn), i != j
     * U = U' + 2J
     * J = J'
     */
    void interaction(int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);
    void exchange(int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);
    void pairHopping(int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap);

    void push(pairIdx_t pairRepIf, T val, MAP<T>* rowMap);
    
protected:
    Basis* pt_Basis{nullptr};
    LATTICE_MODEL fmodel{LATTICE_MODEL::HUBBARD};
    int nuSign;
    bool useTrans{true};
    bool usePG{true};
};

template<typename T>
class SpinOperator{
/*
    idx:0 --> spin:s
    idx++ ---> s--
    For tJ model:
        Spin-1/2 ops for tJ model. Operate on projected hilbert space without double occupancy.(docc)
        In order to not consider fermion sign for spin operators, the basis state is aranged as [su0,sd0,su1,sd1,...]
        The basis state is stored as [su0,su1,...], [sd0,sd1,...]
*/
public:
    SpinOperator( ) { }
    SpinOperator(Basis* pt_Ba, bool commuteWithTrans = false, bool commuteWithPG = false);
    ~SpinOperator( ) { }

    void push(idx_t repIf, T val, MAP<T>* rowMap);

    /*
        VecI Reps
    */
    double getSz(int siteI, VecI& initVec) const { return szMat.at(initVec.at(siteI)); }
    
    void szsz(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap);
    
    void spsm(int siteI, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap);
    
    void smsp(int siteI, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap);
    
    void spsm(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap);
    
    void smsp(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap);
    /*
        Binary Reps For Spindim=2
    */
    double getSz(int siteI, idx_t repI) const;

    void szsz(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);

    void szsznn(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);
    
    void spsm(int siteI, T factor, idx_t repI, MAP<T>* rowMap);
    
    void smsp(int siteI, T factor, idx_t repI, MAP<T>* rowMap);
   
    void spsm(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);
 
    void smsp(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);

    void szspsm(int siteI, int siteJ, int siteK, T factor, idx_t repI, MAP<T>* rowMap);

    void chiral(int siteI, int siteJ, int siteK, T factor, idx_t repI, MAP<T>* rowMap);

protected:
    Basis* pt_Basis{nullptr};
    LATTICE_MODEL smodel{LATTICE_MODEL::HEISENBERG};
    int spinDim{2};
    std::vector<double> szMat, spMat, smMat;
    bool useTrans{true};
    bool usePG{true};
};

/********************
 * FERMION OPERATOR *
 ********************/

template <typename T>
FermionOperator<T>::FermionOperator(Basis* pt_Ba, bool Trans, bool PG):pt_Basis(pt_Ba),fmodel(pt_Ba->getModel()), useTrans(!Trans), usePG(!PG) { }

template <typename T>
bool FermionOperator<T>::cp(SPIN spin, int siteI, pairIdx_t& pairRepI, int &sign){
    switch(fmodel){
        case LATTICE_MODEL::HUBBARD: {
            if (spin == SPIN::UP) {
                if (!bitTest(pairRepI.first,siteI)) {
                    bitFlip(pairRepI.first,siteI);
                    int counter = 0;
                    for (int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.first,i)) ++counter;
                    }
                    int tmp = (counter%2==0)?1:-1;
                    sign *= tmp;
                    return true;
                }
            } else {
                if (!bitTest(pairRepI.second,siteI)) {
                    auto N = pt_Basis->getOrbNum();
                    int counter = 0; 
                    for(int i=0; i<N; ++i) {
                        if (bitTest(pairRepI.first,i)) ++counter;
                    }
                    bitFlip(pairRepI.second,siteI);
                    for (int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.second,i)) ++counter;
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            }
            break;
        }
        case LATTICE_MODEL::tJ: {
            if (spin == SPIN::UP) {
                if (!bitTest(pairRepI.first, siteI) and !bitTest(pairRepI.second, siteI)) {
                    bitFlip(pairRepI.first, siteI);
                    int counter = 0;
                    for (int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.first, i)) {
                            ++counter;
                        } else if (bitTest(pairRepI.second, i)) {
                            ++counter;
                        }
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            } else {
                if (!bitTest(pairRepI.second, siteI) and !bitTest(pairRepI.first, siteI)) {
                    bitFlip(pairRepI.second, siteI);
                    int counter = 0; 
                    for(int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.second, i)) {
                            ++counter;
                        } else if (bitTest(pairRepI.first, i)) {
                            ++counter;
                        }
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            }
            break;
        }
        default:
            std::cout<<"cp not defined for "<<fmodel<<"\n";
            exit(1);
    }
    return false;
}

template <typename T>
bool FermionOperator<T>::cm(SPIN spin, int siteI, pairIdx_t& pairRepI, int &sign){
    switch (fmodel) {
        case LATTICE_MODEL::HUBBARD: {
            if (spin == SPIN::UP) {
                if (bitTest(pairRepI.first, siteI)) {
                    bitFlip(pairRepI.first, siteI);
                    int counter = 0;
                    for(int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.first, i)) ++counter;
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            } else {
                if (bitTest(pairRepI.second, siteI)) {
                    auto N = pt_Basis->getOrbNum();
                    int counter = 0; 
                    for (int i = 0; i <N ; ++i) {
                        if (bitTest(pairRepI.first, i)) counter++;
                    }
                    bitFlip(pairRepI.second, siteI);
                    for(int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.second, i)) ++counter;
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            }
            break;
        }
        case LATTICE_MODEL::tJ: {
            if (spin == SPIN::UP) {
                if (bitTest(pairRepI.first, siteI)) {
                    bitFlip(pairRepI.first, siteI);
                    int counter = 0;
                    for (int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.first, i)) {
                            ++counter;
                        } else if (bitTest(pairRepI.second, i)) {
                            ++counter;
                        }
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            } else {
                if (bitTest(pairRepI.second, siteI)) {
                    bitFlip(pairRepI.second, siteI);
                    int counter = 0; 
                    for(int i = 0; i < siteI; ++i) {
                        if (bitTest(pairRepI.second, i)) {
                            ++counter;
                        } else if (bitTest(pairRepI.first, i)) {
                            ++counter;
                        }
                    }
                    int tmp = (counter % 2 == 0) ? 1 : -1;
                    sign *= tmp;
                    return true;
                }
            } 
            break;
        }
        default:
            std::cout<<"cm not defined for "<<fmodel<<"\n";
            exit(1);
    }
    return false;
}

template <typename T>
bool FermionOperator<T>::cpcm(SPIN spin, int siteI, int siteJ, pairIdx_t& pairRepI, int& sign){
    idx_t& repI = (spin == SPIN::UP) ? pairRepI.first : pairRepI.second;
    if (siteI==siteJ){
        if (bitTest(repI, siteJ)){
            sign = 1;
            return true;
        }
    } else {
        switch(fmodel){
            case LATTICE_MODEL::HUBBARD: {
                if (bitTest(repI, siteJ) && (!bitTest(repI, siteI))){
                    bitFlip(repI, siteI);
                    bitFlip(repI, siteJ);
                    int counter = 0;
                    if (siteI < siteJ) {
                        for (int i = siteI + 1; i < siteJ; ++i) {
                            if (bitTest(repI,i)) ++counter;
                        }
                    } else {
                        for (int i = siteJ + 1; i < siteI; ++i) {
                            if (bitTest(repI,i)) ++counter;
                        }
                    }
                    sign = (counter % 2 == 0) ? 1 : -1;
                    return true;
                }
                break;
            }

            // no double occupancy
            case LATTICE_MODEL::tJ: {
                idx_t repIp = (spin == SPIN::DOWN) ? pairRepI.first : pairRepI.second;
                if (bitTest(repI, siteJ) && (!bitTest(repI, siteI)) && (!bitTest(repIp, siteI))) {
                    bitFlip(repI,siteI);
                    bitFlip(repI,siteJ);
                    int counter = 0;
                    if (siteI < siteJ) {
                        for (int i = siteI + 1; i < siteJ; ++i) {
                            if (bitTest(repI,i)) {
                                ++counter; 
                            } else if (bitTest(repIp,i)) {
                                ++counter;
                            }
                        }
                    }else{
                        for (int i = siteJ + 1; i < siteI; ++i) {
                            if (bitTest(repI,i)) {
                                ++counter;
                            } else if (bitTest(repIp,i)) {
                                ++counter;
                            }
                        }
                    }
                    sign = (counter % 2 == 0) ? 1 : -1;
                    return true;
                }
                break;
            }
            default:break;
        }
    }
    return false;
}

template <typename T>
void FermionOperator<T>::cp(SPIN spin, int siteI, T factor, pairIdx_t pairRepI, MAP<T>* rowMap) {
    int sign = 1;
    idx_t colidx;
    if (cp(spin, siteI, pairRepI, sign)) {
        if (pt_Basis->search(pairRepI, colidx)) {
            T val = factor * (double)sign;
            val /= pt_Basis->getNorm(colidx);
            MapPush(rowMap,colidx,val);
        }
    }
}

template <typename T>
void FermionOperator<T>::cm(SPIN spin, int siteI, T factor, pairIdx_t pairRepI, MAP<T>* rowMap){
    int sign = 1;
    idx_t colidx;
    if (cm(spin, siteI, pairRepI, sign)) {
        if (pt_Basis->search(pairRepI, colidx)) {
            T val = factor * (double)sign;
            val /= pt_Basis->getNorm(colidx);
            MapPush(rowMap,colidx,val);
        }
    }
}

template <typename T>
void FermionOperator<T>::cpcm(SPIN spin, int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap){
    int sign = 1;
    if (cpcm(spin, siteI, siteJ, pairRepI, sign)) {
        push(pairRepI, factor * (double)sign, rowMap);
    }
}

/**
 * @brief Interband exchange interaction: J * d_dag_j_si * d_dag_i_sj * d_j_sj * d_i_si, i != j
 * 
 * @tparam T 
 * @param siteI 
 * @param siteJ 
 * @param factor 
 * @param repI 
 * @param rowMap 
 */
template <typename T>
void FermionOperator<T>::exchange(int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap){
    if (siteI == siteJ) return;
    std::vector<SPIN> spins{SPIN::UP,SPIN::DOWN};
    for(auto spinI:spins){
        for(auto spinJ:spins){
            auto pairRepIf = pairRepI;
            int sign = 1;
            if(cm(spinI, siteI, pairRepIf, sign)){
                if(cm(spinJ, siteJ, pairRepIf, sign)){
                    if(cp(spinJ, siteI, pairRepIf, sign)){
                        if(cp(spinI, siteJ, pairRepIf, sign)){
                            push(pairRepIf, factor * (double)sign, rowMap); 
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Pair hopping: J' * d_dag_i_up * d_dag_i_dn * d_j_dn * d_j_up, i != j
 *  https://doi.org/10.1016/S0370-1573(00)00121-6
 * 
 * @tparam T 
 * @param siteI 
 * @param siteJ 
 * @param factor 
 * @param repI 
 * @param rowMap 
 */
template <typename T>
void FermionOperator<T>::pairHopping(int siteI, int siteJ, T factor, pairIdx_t pairRepI, MAP<T>* rowMap){
    if (siteI == siteJ) return;
    auto pairRepIf = pairRepI;
    int sign = 1;
    if(cm(SPIN::UP, siteI, pairRepIf, sign)){
        if(cm(SPIN::DOWN, siteI, pairRepIf, sign)){
            if(cp(SPIN::DOWN, siteJ, pairRepIf, sign)){
                if(cp(SPIN::UP, siteJ, pairRepIf, sign)){
                    push(pairRepIf, factor * (double)sign, rowMap);
                }
            }
        }
    }

    pairRepIf = pairRepI;
    sign = 1;
    if(cm(SPIN::UP, siteJ, pairRepIf, sign)){
        if(cm(SPIN::DOWN, siteJ, pairRepIf, sign)){
            if(cp(SPIN::DOWN, siteI, pairRepIf, sign)){
                if(cp(SPIN::UP, siteI, pairRepIf, sign)){
                    push(pairRepIf, factor * (double)sign, rowMap);
                }
            }
        }
    }
}

template <typename T>
void FermionOperator<T>::push(pairIdx_t pairRepIf, T val, MAP<T>* rowMap) {
    #ifdef DISTRIBUTED_BASIS
    if(pt_Basis->isfMin(pairRepIf.first)) MapPush(rowMap,pt_Basis->getRepI(pairRepIf),val);
    #else
    idx_t colidx;
    T fac;
    if (pt_Basis->search(pairRepIf, colidx, fac, useTrans, usePG)){
        val /= pt_Basis->getNorm(colidx);
        val *= fac;
        MapPush(rowMap,colidx,val);
    }
    #endif // DISTRIBUTED_BASIS 
}

#ifdef DISTRIBUTED_BASIS

template <typename T>
void FermionOperator<T>::diag(idx_t rowID, T factor, MAP<T>* rowMap){
    idx_t repI = pt_Basis->getRepI(rowID);
    MapPush(rowMap, repI, factor);
}

#else

template <typename T>
void FermionOperator<T>::diag(idx_t rowID, T factor, MAP<T>* rowMap){
        MapPush(rowMap, rowID, factor);
}

#endif // DISTRIBUTED_BASIS


/*****************
 * SPIN OPERATOR *
 *****************/
template<typename T>
SpinOperator<T>::SpinOperator(Basis* pt_Ba, bool Trans, bool PG):\
pt_Basis(pt_Ba), smodel(pt_Ba->getModel()), spinDim(pt_Ba->getSiteDim()), useTrans(!Trans), usePG(!PG) {
    double s = (double)(spinDim - 1)/2.0;
    double m = s;
    for (int i = 0; i < spinDim; ++i){
        szMat.push_back(m);
        spMat.push_back(std::sqrt(s*(s+1.0)-m*(m+1.0)));
        smMat.push_back(std::sqrt(s*(s+1.0)-m*(m-1.0)));
        m -= 1.0;
    }
}

template<typename T>
void SpinOperator<T>::push(idx_t repIf, T val, MAP<T>* rowMap){
    #ifdef DISTRIBUTED_BASIS
        MapPush(rowMap,repIf,val);
    #else
        idx_t colID;
        cdouble fac;
        //!Test
        if (pt_Basis->search(repIf, colID, fac, false, false)){
            double finalNorm = pt_Basis->getNorm(colID);
            val /= finalNorm;
            val *= fac;
            MapPush(rowMap,colID,val);
        }
    #endif
}

template<typename T>
void SpinOperator<T>::szsz(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap){
    #ifdef DISTRIBUTED_BASIS
        T dval = factor * szMat.at(initVec[siteI]) * szMat.at(initVec[siteJ]);
        MapPush(rowMap,initInd,dval);
    #else
        idx_t colID;
        if (pt_Basis->search(initInd,colID)){
            T dval = factor * szMat.at(initVec[siteI]) * szMat.at(initVec[siteJ]);
            double finalNorm = pt_Basis->getNorm(colID);
            dval /= finalNorm;
            MapPush(rowMap,colID,dval);
        }
    #endif
}

template<typename T>
void SpinOperator<T>::spsm(int siteI, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap){
    #ifdef DISTRIBUTED_BASIS
        if (initVec[siteI] < (spinDim-1)){
            T val = factor * spMat[initVec[siteI]+1] * smMat[initVec[siteI]];
            MapPush(rowMap,initInd,val);
        }
    #else
        idx_t colID;
        if (pt_Basis->search(initInd, colID)){
            if (initVec[siteI] < (spinDim-1)){
                T val = factor * spMat[initVec[siteI]+1] * smMat[initVec[siteI]];
                double finalNorm = pt_Basis->getNorm(colID);
                val /= finalNorm;
                MapPush(rowMap,colID,val);
            }
        }
    #endif
}

template<typename T>
void SpinOperator<T>::smsp(int siteI, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap){
    #ifdef DISTRIBUTED_BASIS
        if (initVec[siteI] > 0){
            T val = factor * smMat[initVec[siteI]-1] * spMat[initVec[siteI]];
            MapPush(rowMap,initInd,val);
        }
    #else
        idx_t colID;
        if (pt_Basis->search(initInd, colID)){
            if (initVec[siteI] > 0){
                T val = factor * smMat[initVec[siteI]-1] * spMat[initVec[siteI]];
                double finalNorm = pt_Basis->getNorm(colID);
                val /= finalNorm;
                MapPush(rowMap,colID,val);
            }
        }
    #endif
}

template<typename T>
void SpinOperator<T>::spsm(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap){
    if (siteI==siteJ){
        spsm(siteI, factor, initInd, initVec, rowMap);
        return;
    }
    if (initVec[siteI] > 0 and initVec[siteJ] < (spinDim-1)){
        idx_t finalInd = initInd - pt_Basis->getmul(siteI) + pt_Basis->getmul(siteJ);
        #ifdef Distributde_Basis
            T val = factor * spMat[initVec[siteI]] * smMat[initVec[siteJ]];
            MapPush(rowMap,initInd,val);
        #else
            idx_t colID;
            if (pt_Basis->search(finalInd, colID)){
                T val = factor * spMat[initVec[siteI]] * smMat[initVec[siteJ]];
                double finalNorm = pt_Basis->getNorm(colID);
                val /= finalNorm;
                MapPush(rowMap,colID,val);
            }
        #endif
    }
}

template<typename T>
void SpinOperator<T>::smsp(int siteI, int siteJ, T factor, idx_t initInd, VecI& initVec, MAP<T>* rowMap){
    if (siteI==siteJ){
        smsp(siteI, factor, initInd, initVec, rowMap);
        return;
    }
    spsm(siteJ, siteI, factor, initInd, initVec,rowMap);
}
/**************
 * BINARY REP *
 **************/
template<typename T>
double SpinOperator<T>::getSz(int siteI, idx_t repI) const {
    switch(smodel){
        case LATTICE_MODEL::HEISENBERG:{
            return szMat.at(1&(repI>>siteI));
            break;
        }
        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI=pt_Basis->getPairRepI(repI);
            if(bitTest(pairRepI.first,siteI)){
                return szMat.at(0);
            }else if(bitTest(pairRepI.second,siteI)){
                return szMat.at(1);
            }else{
                return 0.0;
            }
            break;
        }
        default:{
            std::cout<<"model not defined for SpinOperator::getSz\n";
            exit(1);
        }
    }   
}

template<typename T>
void SpinOperator<T>::szsz(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    #ifdef DISTRIBUTED_BASIS
        T dval = factor * getSz(siteI,repI) * getSz(siteJ,repI);
        MapPush(rowMap,repI,dval);
    #else
        T dval = factor * getSz(siteI,repI) * getSz(siteJ,repI);
        push(repI, dval, rowMap);
    #endif
}

template<typename T>
void SpinOperator<T>::szsznn(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    // for tJ model, szi*szj -1/4*ni*nj
    assert_msg(smodel==LATTICE_MODEL::tJ,"SpinOperator::szsznn only defined for tJ model");
    pairIdx_t pairRepI = pt_Basis->getPairRepI(repI);
    if((bitTest(pairRepI.first,siteI) && bitTest(pairRepI.second,siteJ)) || (bitTest(pairRepI.first,siteJ) && bitTest(pairRepI.second,siteI))){
        #ifdef DISTRIBUTED_BASIS
            T dval = -0.5*factor;
            MapPush(rowMap,repI,dval);
        #else
            T dval = -0.5*factor;
            push(repI, dval, rowMap);
        #endif
    }
}

template<typename T>
void SpinOperator<T>::spsm(int siteI, T factor, idx_t repI, MAP<T>* rowMap){
    bool condition;
    switch(smodel){
        case LATTICE_MODEL::HEISENBERG:{
            condition = !bitTest(repI,siteI);
            break;
        }
        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI = pt_Basis->getPairRepI(repI);
            condition = bitTest(pairRepI.first,siteI);
            break;
        }
        default:{
            std::cout<<"model not defined for SpinOperator::spsm\n";
            exit(1);
        }
    }
    if (condition){
        #ifdef DISTRIBUTED_BASIS
            T val = factor * spMat[1] * smMat[0];
            MapPush(rowMap,repI,val);
        #else
            T val = factor * spMat[1] * smMat[0];
            push(repI, val, rowMap);
        #endif
    }
}

template<typename T>
void SpinOperator<T>::smsp(int siteI, T factor, idx_t repI, MAP<T>* rowMap){
    bool condition;
    switch(smodel){
        case LATTICE_MODEL::HEISENBERG:{
            condition = bitTest(repI,siteI);
            break;
        }
        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI = pt_Basis->getPairRepI(repI);
            condition = bitTest(pairRepI.second,siteI);
            break;
        }
        default:{
            std::cout<<"model not defined for SpinOperator::smsp\n";
            exit(1);
        }
    }
    if (condition){
        #ifdef DISTRIBUTED_BASIS
            T val = factor * spMat[1] * smMat[0];
            MapPush(rowMap,repI,val);
        #else
            T val = factor * spMat[1] * smMat[0];
            push(repI, val, rowMap);
        #endif
    }
}

template<typename T>
void SpinOperator<T>::spsm(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    if (siteI==siteJ){
        spsm(siteI, factor, repI, rowMap);
        return;
    }
    switch(smodel){
        case LATTICE_MODEL::HEISENBERG:{
            if (bitTest(repI,siteI) && (!bitTest(repI,siteJ))){
                bitFlip(repI,siteI);
                bitFlip(repI,siteJ);
                T val = factor * spMat[1] * smMat[0];
                #ifdef DISTRIBUTED_BASIS
                    MapPush(rowMap,repI,val);
                #else
                    push(repI, val, rowMap);
                #endif
            }
            break;
        }
        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI = pt_Basis->getPairRepI(repI);
            if (bitTest(pairRepI.first,siteJ) && bitTest(pairRepI.second,siteI)){
                bitFlip(pairRepI.first,siteI);
                bitFlip(pairRepI.first,siteJ);
                bitFlip(pairRepI.second,siteI);
                bitFlip(pairRepI.second,siteJ);
                repI = pt_Basis->getRepI(pairRepI);
                T val = factor * spMat[1] * smMat[0];
                #ifdef DISTRIBUTED_BASIS
                    MapPush(rowMap,repI,val);
                #else
                push(repI, val, rowMap);
                #endif
            }
            break;
        }
        default:{
            std::cout<<"model not defined for SpinOperator::spsm\n";
            exit(1);
        }
    }
}

template<typename T>
void SpinOperator<T>::smsp(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    if (siteI==siteJ){
        smsp(siteI, factor, repI,rowMap);
        return;
    }
    spsm(siteJ, siteI, factor, repI,rowMap);
}

template<typename T>
void SpinOperator<T>::szspsm(int siteI, int siteJ, int siteK, T factor, idx_t repI, MAP<T>* rowMap){
    /*
        Jk/norm(repI)/norm(repfI) * I/2 * sz(i) * [sp(j)sm(k) - sm(j)sp(k)]
        factor = Jk/norm(repI)
    */
    // assert_msg(siteI!=siteJ && siteJ!=siteK, "chiral term szspsm(i,j,k) only defined for three different sites (i,j,k)!");
    switch (smodel){
        case LATTICE_MODEL::HEISENBERG:{
            bool cond = false;
            cdouble sign;
            if (bitTest(repI,siteJ) && (!bitTest(repI,siteK))){
                cond = true;
                sign = CPLX_I/2.0;
            }
            else if ((!bitTest(repI,siteJ)) && bitTest(repI,siteK)){
                cond = true;
                sign = -CPLX_I/2.0;
            }
            if (cond){
                bitFlip(repI,siteJ);
                bitFlip(repI,siteK);
                #ifdef DISTRIBUTED_BASIS
                    T val = factor * (double)sign * getSz(siteI, repI) * spMat[1] * smMat[0];
                    MapPush(rowMap,repI,val);
                #else
                    T val = sign * factor * getSz(siteI, repI) * spMat[1] * smMat[0];
                    push(repI, val, rowMap);
                #endif
            }
            break;
        }
        default:{
            std::cout<<"model not defined for SpinOperator::szspsm\n";
            exit(1);
        }
    }
}

template<typename T>
void SpinOperator<T>::chiral(int siteI, int siteJ, int siteK, T factor, idx_t repI, MAP<T>* rowMap){
    /*
        Jk/norm(repI)/norm(repfI) * si * (sj x sk)
        factor = Jk/norm(repI)
    */
    assert_msg(siteI!=siteJ && siteJ!=siteK, "chiral term chiral(i,j,k) only defined for three different lattice sites (i,j,k)!");
    szspsm(siteI, siteJ, siteK, factor, repI, rowMap);
    szspsm(siteJ, siteK, siteI, factor, repI, rowMap);
    szspsm(siteK, siteI, siteJ, factor, repI, rowMap);
}