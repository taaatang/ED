#ifndef __OPERATORBASICS_H__
#define __OPERATORBASICS_H__

#include "Basis/Basis.hpp"

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
public:
    Link( ) { }
    Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, bool is_Const_ = true, bool is_Ordered_ = false);
    Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, const std::vector<VecD>& vecs, bool is_Const_ = true, bool is_Ordered_ = false);
    ~Link( ) { };
    
    LINK_TYPE getLinkType( ) const { return link_type; }
    T getVal( ) const { return val; }
    int getlinkid( ) const { return linkid; }
    int getmatid( ) const { return matid; }
    int getLinkNum( ) const { return LinkList.size(); }
    int getvecid(int bondid) const { return LinkVecIdList.at(bondid); }
    VecD getvec(int vecid) const { return LinkVecs.at(vecid); }
    // rep orbital of the link
    ORBITAL getRepOrb( ) const { return orbList.at(0); }
    std::vector<ORBITAL> getOrbs( ) const { return orbList; }
    // number of orbs in one link
    int getLinkSize( ) const { return orbList.size() - 1; }
    // size of LinkVecs
    int getLinkVecNum( ) const { return LinkVecs.size(); }

    // is this a const link or time dependent link
    bool isConst( ) const { return is_Const; }
    bool isOrdered( ) const { return is_Ordered; }
    void setConst(bool cond) { is_Const = cond; }
    void setOrdered(bool cond) { is_Ordered = cond; }

    // begin/end iterator for LinkList
    std::vector<VecI>::iterator begin( ) { return LinkList.begin(); }
    std::vector<VecI>::iterator end( ) { return LinkList.end(); }
    const std::vector<VecI>& bond( ) const { return LinkList; }

    void setVal(T val_) { val = val_; }
    void setid(int linkid_, int matid_ = 0);
    void setmatid(int matid_ = 0) { matid = matid_;}
    // set time evolve function for val
    void setTimeFunc(T (*timeFunc_)(int));
    void setVal( );

    // add one link (contains orbids in the same link)
    void push_back(std::vector<int> orbidList) { LinkList.push_back(orbidList); }
    // add a orb's relative coord to the rep orb
    Link<T>& addLinkVec(std::vector<double> vec);
    //generate all the links->linkList
    void genLinkMaps(Geometry* pt_lattice);
    void print(bool brief = true, std::ostream& os = std::cout) const;

private:
    LINK_TYPE link_type{LINK_TYPE::HOPPING_T};
    T val{0.0};
    // same group ids will be contained in the same sparse matrix
    int linkid{0}, matid{0};
    bool is_Const{true}, is_Ordered{false};
    T (*timeFunc)(int);
    bool is_timeFunc_set{false};
    int timeStep{0};
    std::vector<ORBITAL> orbList;
    std::vector<VecD> LinkVecs;
    std::vector<VecI> LinkList;
    // LinkVecIdList[i] is the id of LinkVec for i-th bond in LinkList
    VecI LinkVecIdList;
};

/*
    *******************
    * Operators Class *
    *******************
*/
template<typename T>
class FermionOperator{
public:
    FermionOperator( ) { }
    FermionOperator(Basis* pt_Ba);
    ~FermionOperator( ) { };
    // c^dagger/c act on siteI of repI
    bool cp(int siteI, idx_t repI, idx_t& repIf, int &sign);
    bool cm(int siteI, idx_t repI, idx_t& repIf, int &sign);
    // c^dagger/c with spin act on siteI of pairRepI 
    bool cp(SPIN spin, int siteI, pairIndex& pairRepI, int &sign);
    bool cm(SPIN spin, int siteI,  pairIndex& pairRepI, int &sign);
    // used for single particle spectra 
    void cp(SPIN spin, int siteI, T factor, pairIndex pairRepI, MAP<T>* rowMap);
    void cm(SPIN spin, int siteI, T factor, pairIndex pairRepI, MAP<T>* rowMap);
    // hopping deals with sign count difference for hubbard and tJ
    bool cpcm(int siteI, int siteJ, idx_t repI, idx_t repIp, idx_t& repIf, int& sign);
    // hopping term
    void cpcm(SPIN spin, int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);
   
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
    void interaction(int siteI, int siteJ, T factor, idx_t rep_I, MAP<T>* rowMap);
    void exchange(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);
    void pairHopping(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap);
protected:
    Basis* pt_Basis{nullptr};
    LATTICE_MODEL fmodel{HUBBARD};
    int nuSign;
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
    SpinOperator(Basis* pt_Ba);
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
    LATTICE_MODEL smodel{HEISENBERG};
    int spinDim{2};
    std::vector<double> szMat, spMat, smMat;
};


/********
 * LINK *
 ********/
template <class T>
Link<T>::Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, bool is_Const_ , bool is_Ordered_ ):\
timeStep(0), is_timeFunc_set(false) {
    link_type = link_type_;
    orbList = orbList_;
    val = val_;
    is_Const = is_Const_;
    is_Ordered = is_Ordered_;
}

template <class T>
Link<T>::Link(LINK_TYPE link_type_, std::vector<ORBITAL> orbList_, T val_, const std::vector<VecD>& vecs, bool is_Const_ , bool is_Ordered_ ):\
timeStep(0), is_timeFunc_set(false) {
    link_type = link_type_;
    orbList = orbList_;
    val = val_;
    is_Const = is_Const_;
    is_Ordered = is_Ordered_;
    for (auto& vec : vecs) {
        addLinkVec(vec);
    }
}

template <class T>
void Link<T>::setid(int linkid_, int matid_ ) { 
    linkid = linkid_; 
    matid = matid_; 
}

template <class T>
void Link<T>::setTimeFunc(T (*timeFunc_)(int)) {
    timeFunc = timeFunc_; 
    is_timeFunc_set = true;
}

template <class T>
void Link<T>::setVal( ) {
    if(!is_Const) {
        if(!is_timeFunc_set){std::cout<<"time evolve function not set for linkid:"<<linkid<<std::endl; exit(1);}
        val = timeFunc(timeStep); 
    }
    timeStep++;
}

template <class T>
Link<T>& Link<T>::addLinkVec(std::vector<double> vec){
    assert(vec.size()==3); 
    LinkVecs.push_back(vec); 
    return *this;
}

template <class T>
void Link<T>::genLinkMaps(Geometry* pt_lattice) {
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
                    if (pt_lattice->coordToOrbid(orbList.at(k-j+1), coordf.data(), orbidf)){
                        assert(pt_lattice->getOrb(orbidf) == orbList.at(k-j+1));
                        tmp.push_back(orbidf);
                    }else{
                        is_bond = false;
                        break;
                    }
                }
                if(is_bond){
                    push_back(tmp);
                    LinkVecIdList.push_back(j);
                }
            }
        }   
    }
}

template <class T>
void Link<T>::print(bool brief, std::ostream& os) const {
    os<<"mat id:"<<matid<<", Link id:"<<linkid<<", value:"<<val<<"\n";
    for (auto orb:orbList) os<<orb<<" ";
    os<<"\n";
    os<<"link vec:\n";
    for (const auto& vec:LinkVecs) os<<vec<<"\n";
    os<<"bond num:"<<LinkList.size()<<"\n";
    if (!brief) {
        for (const auto& bond : LinkList) os<<bond<<"\n";
    }
}

/********************
 * FERMION OPERATOR *
 ********************/

template <typename T>
FermionOperator<T>::FermionOperator(Basis* pt_Ba):pt_Basis(pt_Ba),fmodel(pt_Ba->getModel()) {

}

template <typename T>
bool FermionOperator<T>::cp(int siteI, idx_t repI, idx_t& repIf, int &sign){
    if(!bitTest(repI,siteI)){
        repIf = repI;
        bitFlip(repIf,siteI);
        int counter = 0;
        for (int i = 0; i < siteI; ++i) {
            if (bitTest(repI,i)) ++counter;
        }
        sign = (counter%2==0)?1:-1;
        return true;
    }
    return false;
}

template <typename T>
bool FermionOperator<T>::cm(int siteI, idx_t repI, idx_t& repIf, int &sign){
    if(bitTest(repI,siteI)){
        repIf = repI;
        bitFlip(repIf,siteI);
        int counter = 0;
        for (int i = 0; i < siteI; ++i) {
            if (bitTest(repI,i)) ++counter;
        }
        sign = (counter%2==0)?1:-1;
        return true;
    }
    return false;    
}

template <typename T>
bool FermionOperator<T>::cp(SPIN spin, int siteI, pairIndex& pairRepI, int &sign){
    switch(fmodel){
        case LATTICE_MODEL::HUBBARD:{
            if(spin == SPIN_UP){
                if(!bitTest(pairRepI.first,siteI)){
                    bitFlip(pairRepI.first,siteI);
                    int counter = 0;
                    for(int i = 0; i < siteI; ++i)if(bitTest(pairRepI.first,i))counter++;
                    int tmp = (counter%2==0)?1:-1;
                    sign *= tmp;
                    return true;
                }
            }else{
                if(!bitTest(pairRepI.second,siteI)){
                    auto N = pt_Basis->getOrbNum();
                    int counter = 0; for(int i=0;i<N;++i)if(bitTest(pairRepI.first,i))counter++;
                    bitFlip(pairRepI.second,siteI);
                    for(int i = 0; i < siteI; ++i)if(bitTest(pairRepI.second,i))counter++;
                    int tmp = (counter%2==0)?1:-1;
                    sign *= tmp;
                    return true;
                }
            }
            break;
        }
        case LATTICE_MODEL::tJ:{
            if(spin == SPIN_UP){
                if(!bitTest(pairRepI.first,siteI) and !bitTest(pairRepI.second,siteI)){
                    bitFlip(pairRepI.first,siteI);
                    int counter = 0;
                    for(int i = 0; i < siteI; ++i) if(bitTest(pairRepI.first,i)) ++counter;
                    int tmp = (counter%2==0)?1:-1;
                    sign *= tmp;
                    return true;
                }
            }else{
                if(!bitTest(pairRepI.second,siteI) and !bitTest(pairRepI.first,siteI)){
                    auto N = pt_Basis->getOrbNum();
                    int counter = 0; for(int i=0;i<N;++i) if(bitTest(pairRepI.first,i)) ++counter;
                    bitFlip(pairRepI.second,siteI);
                    for(int i = 0; i < siteI; ++i)if(bitTest(pairRepI.second,i))counter++;
                    int tmp = (counter%2==0)?1:-1;
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
bool FermionOperator<T>::cm(SPIN spin, int siteI, pairIndex& pairRepI, int &sign){
    if(spin == SPIN_UP){
        if(bitTest(pairRepI.first,siteI)){
            bitFlip(pairRepI.first,siteI);
            int counter = 0;
            for(int i = 0; i < siteI; ++i)if(bitTest(pairRepI.first,i))counter++;
            int tmp = (counter%2==0)?1:-1;
            sign *= tmp;
            return true;
        }
    }else{
        if(bitTest(pairRepI.second,siteI)){
            auto N = pt_Basis->getOrbNum();
            int counter = 0; for(int i=0;i<N;++i)if(bitTest(pairRepI.first,i))counter++;
            bitFlip(pairRepI.second,siteI);
            for(int i = 0; i < siteI; ++i)if(bitTest(pairRepI.second,i))counter++;
            int tmp = (counter%2==0)?1:-1;
            sign *= tmp;
            return true;
        }
    }
    return false;
}

template <typename T>
void FermionOperator<T>::cp(SPIN spin, int siteI, T factor, pairIndex pairRepI, MAP<T>* rowMap) {
    pairIndex pairRepIf = pairRepI;
    int sign;
    idx_t colidx;
    switch (spin)
    {
        case SPIN_UP:
            if (cp(siteI, pairRepI.first, pairRepIf.first, sign)) {
                if (pt_Basis->search(pairRepIf,colidx)) {
                    T val = factor * (double)sign;
                    val /= pt_Basis->getNorm(colidx);
                    MapPush(rowMap,colidx,val);
                }
            }
            break;
        case SPIN_DOWN:
            if(cp(siteI, pairRepI.second, pairRepIf.second, sign)){
                if(pt_Basis->search(pairRepIf,colidx)){
                    auto N = pt_Basis->getOrbNum();
                    int count = 0; 
                    for (int i = 0; i < N; ++i) {
                        if (bitTest(pairRepI.first,i)) count++;
                    }
                    if (count % 2 == 1) sign *= -1;
                    T val = factor * (double)sign;
                    val /= pt_Basis->getNorm(colidx);
                    MapPush(rowMap,colidx,val);
                }
            }
            break;
        default:
            break;
    }
}

template <typename T>
void FermionOperator<T>::cm(SPIN spin, int siteI, T factor, pairIndex pairRepI, MAP<T>* rowMap){
    pairIndex pairRepIf = pairRepI;
    int sign;
    idx_t colidx;
    switch (spin)
    {
        case SPIN_UP:
            if(cm(siteI, pairRepI.first, pairRepIf.first, sign)){
                if (pt_Basis->search(pairRepIf,colidx)) {
                    T val = factor * (double)sign;
                    val /= pt_Basis->getNorm(colidx);
                    MapPush(rowMap,colidx,val);
                }
            }
            break;
        case SPIN_DOWN:
            if (cm(siteI, pairRepI.second, pairRepIf.second, sign)) {
                if (pt_Basis->search(pairRepIf,colidx)) {
                    auto N = pt_Basis->getOrbNum();
                    int count = 0; 
                    for (int i = 0; i < N; ++i) {
                        if (bitTest(pairRepI.first,i)) count++;
                    }
                    if (count % 2 == 1) sign *= -1;
                    T val = factor * (double)sign;
                    val /= pt_Basis->getNorm(colidx);
                    MapPush(rowMap,colidx,val);
                }
            }
            break;
        default:
            break;
    }
}

template <typename T>
bool FermionOperator<T>::cpcm(int siteI, int siteJ, idx_t repI, idx_t repIp, idx_t& repIf, int& sign){
    repIf = repI;
    if (siteI==siteJ){
        if (bitTest(repI,siteJ)){
            sign = 1;
            return true;
        }
    }else{
        switch(fmodel){
            case LATTICE_MODEL::HUBBARD:{
                if(bitTest(repI,siteJ) && (!bitTest(repI,siteI))){
                    bitFlip(repIf,siteI);
                    bitFlip(repIf,siteJ);
                    int counter = 0;
                    if (siteI < siteJ) {
                        for (int i = siteI + 1; i < siteJ; ++i) if (bitTest(repI,i)) counter++;
                    }else{
                        for (int i = siteJ + 1; i < siteI; ++i) if (bitTest(repI,i)) counter++;
                    }
                    sign = (counter%2==0)?1:-1;
                    return true;
                }
                break;
            }

            // no double occupancy
            case LATTICE_MODEL::tJ:{
                if(bitTest(repI,siteJ) && (!bitTest(repI,siteI)) && (!bitTest(repIp,siteI))){
                    bitFlip(repIf,siteI);
                    bitFlip(repIf,siteJ);
                    int counter = 0;
                    if (siteI < siteJ) {
                        for (int i = siteI + 1; i < siteJ; ++i) {if (bitTest(repI,i)) counter++; else if (bitTest(repIp,i)) counter++;}
                    }else{
                        for (int i = siteJ + 1; i < siteI; ++i) {if (bitTest(repI,i)) counter++; else if (bitTest(repIp,i)) counter++;}
                    }
                    sign = (counter%2==0)?1:-1;
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
void FermionOperator<T>::cpcm(SPIN spin, int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    pairIndex pairRepI=pt_Basis->getPairRepI(repI);
    pairIndex pairRepIf = pairRepI;
    int sign;
    idx_t colidx;
    switch(spin){
        case SPIN_UP:{
            if (cpcm(siteI, siteJ, pairRepI.first, pairRepI.second, pairRepIf.first, sign)){
                #ifdef DISTRIBUTED_BASIS
                if(pt_Basis->isfMin(pairRepIf.first)) MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                #else
                if (pt_Basis->search(pairRepIf, colidx)){
                    T val = factor * (double)sign;
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
                    T val = factor * (double)sign;
                    double finalNorm = pt_Basis->getNorm(colidx);
                    val /= finalNorm;
                    MapPush(rowMap,colidx,val);
                }
                #endif
            }
            break;
        }
        default:
            break;
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
void FermionOperator<T>::exchange(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    if (siteI == siteJ) return;
    auto pairRepI = pt_Basis->getPairRepI(repI);
    idx_t colidx;
    std::vector<SPIN> spins{SPIN_UP,SPIN_DOWN};
    for(auto spinI:spins){
        for(auto spinJ:spins){
            auto pairRepIf = pairRepI;
            int sign = 1;
            if(cm(spinI,siteI,pairRepIf,sign)){
                if(cm(spinJ,siteJ,pairRepIf,sign)){
                    if(cp(spinJ,siteI,pairRepIf,sign)){
                        if(cp(spinI,siteJ,pairRepIf,sign)){
                            #ifdef DISTRIBUTED_BASIS
                            MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                            #else
                            if (pt_Basis->search(pairRepIf, colidx)){
                                T val = factor * (double)sign;
                                double finalNorm = pt_Basis->getNorm(colidx);
                                val /= finalNorm;
                                MapPush(rowMap,colidx,val);
                            }
                            #endif 
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
void FermionOperator<T>::pairHopping(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    if (siteI == siteJ) return;
    auto pairRepIf = pt_Basis->getPairRepI(repI);
    idx_t colidx;
    int sign = 1;
    if(cm(SPIN_UP,siteI,pairRepIf,sign)){
        if(cm(SPIN_DOWN,siteI,pairRepIf,sign)){
            if(cp(SPIN_DOWN,siteJ,pairRepIf,sign)){
                if(cp(SPIN_UP,siteJ,pairRepIf,sign)){
                    #ifdef DISTRIBUTED_BASIS
                    MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                    #else
                    if (pt_Basis->search(pairRepIf, colidx)){
                        T val = factor * (double)sign;
                        double finalNorm = pt_Basis->getNorm(colidx);
                        val /= finalNorm;
                        MapPush(rowMap,colidx,val);
                    }
                    #endif 
                }
            }
        }
    }

    sign = 1;
    if(cm(SPIN_UP,siteJ,pairRepIf,sign)){
        if(cm(SPIN_DOWN,siteJ,pairRepIf,sign)){
            if(cp(SPIN_DOWN,siteI,pairRepIf,sign)){
                if(cp(SPIN_UP,siteI,pairRepIf,sign)){
                    #ifdef DISTRIBUTED_BASIS
                    MapPush(rowMap,pt_Basis->getRepI(pairRepIf),factor*sign);
                    #else
                    if (pt_Basis->search(pairRepIf, colidx)){
                        T val = factor * (double)sign;
                        double finalNorm = pt_Basis->getNorm(colidx);
                        val /= finalNorm;
                        MapPush(rowMap,colidx,val);
                    }
                    #endif 
                }
            }
        }
    }
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
SpinOperator<T>::SpinOperator(Basis* pt_Ba):\
pt_Basis(pt_Ba),smodel(pt_Ba->getModel()),spinDim(pt_Ba->getSiteDim()) {
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
        if (pt_Basis->search(repIf,colID)){
            double finalNorm = pt_Basis->getNorm(colID);
            val /= finalNorm;
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
            pairIndex pairRepI=pt_Basis->getPairRepI(repI);
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
        idx_t colID;
        if (pt_Basis->search(repI,colID)){
            T dval = factor * getSz(siteI,repI) * getSz(siteJ,repI);
            double finalNorm = pt_Basis->getNorm(colID);
            dval /= finalNorm;
            MapPush(rowMap,colID,dval);
        }
    #endif
}

template<typename T>
void SpinOperator<T>::szsznn(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    // for tJ model, szi*szj -1/4*ni*nj
    assert_msg(smodel==LATTICE_MODEL::tJ,"SpinOperator::szsznn only defined for tJ model");
    pairIndex pairRepI = pt_Basis->getPairRepI(repI);
    if((bitTest(pairRepI.first,siteI) && bitTest(pairRepI.second,siteJ)) || (bitTest(pairRepI.first,siteJ) && bitTest(pairRepI.second,siteI))){
        #ifdef DISTRIBUTED_BASIS
            T dval = -0.5*factor;
            MapPush(rowMap,repI,dval);
        #else
            idx_t colID;
            if (pt_Basis->search(repI,colID)){
                T dval = -0.5*factor;
                double finalNorm = pt_Basis->getNorm(colID);
                dval /= finalNorm;
                MapPush(rowMap,colID,dval);
            }
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
            pairIndex pairRepI = pt_Basis->getPairRepI(repI);
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
            idx_t colID;
            if (pt_Basis->search(repI, colID)){
                T val = factor * spMat[1] * smMat[0];
                double finalNorm = pt_Basis->getNorm(colID);
                val /= finalNorm;
                MapPush(rowMap,colID,val);
            }
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
            pairIndex pairRepI = pt_Basis->getPairRepI(repI);
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
            idx_t colID;
            if (pt_Basis->search(repI, colID)){
                T val = factor * spMat[1] * smMat[0];
                double finalNorm = pt_Basis->getNorm(colID);
                val /= finalNorm;
                MapPush(rowMap,colID,val);
            }
        #endif
    }
}

template<typename T>
void SpinOperator<T>::spsm(int siteI, int siteJ, T factor, idx_t repI, MAP<T>* rowMap){
    idx_t colID;
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
                    idx_t colID;
                    if (pt_Basis->search(repI, colID)){
                        double finalNorm = pt_Basis->getNorm(colID);
                        val /= finalNorm;
                        MapPush(rowMap,colID,val);
                    }
                #endif
            }
            break;
        }
        case LATTICE_MODEL::tJ:{
            pairIndex pairRepI = pt_Basis->getPairRepI(repI);
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
                idx_t colID;
                if (pt_Basis->search(repI, colID)){
                    double finalNorm = pt_Basis->getNorm(colID);
                    val /= finalNorm;
                    MapPush(rowMap,colID,val);
                }
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
                    idx_t colID;
                    if (pt_Basis->search(repI, colID)){
                        T val = sign * factor * getSz(siteI, repI) * spMat[1] * smMat[0];
                        double finalNorm = pt_Basis->getNorm(colID);
                        val /= finalNorm;
                        MapPush(rowMap,colID,val);
                    }
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

#endif // __OPERATORBASICS_H__