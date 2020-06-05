//
//  Basis.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Basis.hpp"

/*
    ***************
    * Basis Class *
    ***************
*/
// initialize
Basis::Basis(LATTICE_MODEL input_model, Geometry *pt_lat, VecI& occList, int kInd, int PGRepInd, int siteD):\
    model(input_model), pt_lattice(pt_lat),kIndex(kInd),PGRepIndex(PGRepInd),siteDim(siteD){
    /*
     siteDim is the Hilbert space dimension of a single site. usually equals to spin dimension: 2*s+1
     occList is of length sitedim, containing the number of sites in the state of corresponding site Hilbert Space dimension.
     HUBBARD MODEL: occList->{#spin-up el, #spin-down el}
     t-J Model: occList->{#holes, #spin-up, #spin-down}
     HEISENBERG: occList->{#spin-up, #spin-down}
    */
    #ifdef BINARY_REP
        assert(siteDim==2);
    #endif
    N = pt_lattice->getOrbNum();
    assert(kIndex<pt_lattice->getSiteNum());
    assert(PGRepIndex<pt_lattice->getPGRepNum());
    switch(model){
        case LATTICE_MODEL::HUBBARD:
            Nocc = occList;
            assert(Nocc.at(0)<=N);
            assert(Nocc.at(1)<=N);
            fDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc[0]));
            sDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc[1]));
            totDim = fDim * sDim;
            break;
        case LATTICE_MODEL::t_J:
            Nocc = occList;
            assert((Nocc.at(0)+Nocc.at(1))<=N);
            fDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc[0]));
            sDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc[1]));
            totDim = fDim * combination<ind_int>((ind_int)(N-Nocc[0]), (ind_int)(Nocc[1]));
            break;
        case LATTICE_MODEL::HEISENBERG:
            Sztot = occList.at(0) - occList.at(1);
            Nocc = occList;;
            #ifdef BINARY_REP
                assert((Nocc.at(0) + Nocc.at(1))==N);
                // calculate totDim
                totDim = 1;
                ind_int n = N;
                for (int i = 0; i < siteDim; i++){
                    totDim *= combination<ind_int>(n, (ind_int)(Nocc[i]));
                    n -= Nocc[i];
                }
            #else
                // initialize vec
                for (ind_int i = 0; i < siteDim; i++){
                    for (int j = 0; j < occList[i]; j++) vec.push_back(i);
                }
            #endif
            break;
        default:break;
    }
    initMinMaxRep();
    // default full hilbert space
    subDim = totDim;
    locDim = subDim;
    // initialize mul
    for (int i = 0; i < N; i++) mul.push_back(pow((ind_int) siteDim, (ind_int) (N - i - 1)));  
    // generate fIndexList and sIndexList for Hubbard and tJ
    gendcmp(); 
}
void Basis::initMinMaxRep() const {
#ifdef BINARY_REP
    assert(siteDim==2);
    if(model==LATTICE_MODEL::HUBBARD or model==LATTICE_MODEL::t_J){
        fminRep=0;fmaxRep=0;
        for(int i=0; i<Nocc.at(0); i++){
            bitSet(fminRep,i);
            bitSet(fmaxRep,N-1-i);
        }
        sminRep=0;smaxRep=0;
        for(int i=0; i<Nocc.at(1); i++){
            bitSet(sminRep,i);
            bitSet(smaxRep,N-1-i);
        }
    }
    if(model==LATTICE_MODEL::HEISENBERG){
        fminRep=0;fmaxRep=0;
        for(int i=0; i<Nocc.at(1); i++){
            bitSet(fminRep,i);
            bitSet(fmaxRep,N-1-i);
        }
    }
#else
    vec.clear();
    vecp.clear();
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::t_J:
            for (ind_int i = 0; i < (N - Nocc); i++) vec.push_back(0);
            for (ind_int i = 0; i < Nocc; i++) vec.push_back(1);
            for (ind_int i = 0; i < (Np - Npocc); i++) vecp.push_back(0);
            for (ind_int i = 0; i < Npocc; i++) vecp.push_back(1);
            break;
        default:break;
    }
#endif
}

void Basis::gendcmp(){
    if(model==LATTICE_MODEL::HUBBARD || model==LATTICE_MODEL::t_J){
        fIndexList.clear(); sIndexList.clear();
        fIndexList.reserve(fDim); sIndexList.reserve(sDim);

        ind_int counter = 0;
        ind_int repI = fminRep;
        while(repI<=fmaxRep){
            fIndexList.push_back(repI);
            fRepIdxHash[repI] = counter;

            repI = nextLexicographicalNumber<ind_int>(repI);
            counter++;
        }
        assert(counter == fDim);

        counter = 0;
        repI = sminRep;
        while(repI<=smaxRep){
            sIndexList.push_back(repI);
            sRepIdxHash[repI] = counter;

            repI = nextLexicographicalNumber<ind_int>(repI);
            counter++;
        }
        assert(counter == sDim);

        for(ind_int idx=0; idx<fDim; idx++){
            ind_int repI = fIndexList.at(idx);
            VecI symmList;
            if(isMin(repI, symmList)){
                fMinRepSymmHash[repI] = symmList;
                fMinRepList.push_back(repI);
            } 
        }
    }
}
// generate Basis for the subspace labeled by kInd
void Basis::gen(){
    initMinMaxRep();
    subDim = 0;
    indexList.clear();
    normList.clear();
    ind_int repI;
    double norm;
    switch(model){
        case HUBBARD:case t_J:{
            if (!(kIndex==-1 and model==LATTICE_MODEL::HUBBARD)){ 
                for(ind_int fidx=0;fidx<fDim;fidx++){
                    if(!isfMin(fIndexList[fidx])) continue;
                    for(ind_int sidx=0;sidx<sDim;sidx++){
                        repI = fidx*sDim+sidx;
                        if(isMinRep(repI,norm)){
                            std::cout<<"test repI is minRepI"<<"\n";
                            indexList.push_back(repI);
                            subDim++;
                            #ifdef KEEP_BASIS_NORM
                                normList.push_back(norm);
                            #endif
                        }
                    }
                }
            }else{
                subDim=totDim;
            }
            break;
        }
        case HEISENBERG:{
        #ifdef BINARY_REP
            ind_int repI = fminRep;
            while(repI<=fmaxRep){
                if (isMinRep(repI, norm)){
                    indexList.push_back(repI);
                    subDim++;
                    #ifdef KEEP_BASIS_NORM
                        if (kIndex!=-1) normList.push_back(norm);
                    #endif
                }
                repI = nextLexicographicalNumber<ind_int>(repI);
            }
        #else
            // do{
            //     ind_int repI = vecToRep(vec);
            //     if (isMinRep(repI, norm)){
            //         indexList.push_back(repI);
            //         subDim++;
            //         #ifdef KEEP_BASIS_NORM
            //             if (kIndex!=-1) normList.push_back(norm);
            //         #endif
            //     }
            // }while (std::next_permutation(vec.begin(), vec.end()));
        #endif
            break;
        }
        default:break;                                                  
    }
    if (kIndex==-1) assert(subDim == totDim);
    locDim = subDim;
}

void Basis::gen(int workerID, int workerNum){
    locDim = 0;
    assert(model==LATTICE_MODEL::HUBBARD or model==LATTICE_MODEL::t_J);
    if(model==LATTICE_MODEL::HUBBARD) assert(kIndex!=-1);
    ind_int repI;
    double norm;
    ind_int size = fMinRepList.size();
    ind_int fidxStart, fidxEnd;
    work_load(size, workerID, workerNum, fidxStart, fidxEnd);
    for(ind_int idx=fidxStart; idx<fidxEnd; idx++){
        ind_int fidx = fRepIdxHash.at(fMinRepList[idx]);
        for(ind_int sidx=0;sidx<sDim;sidx++){
            repI = fidx*sDim+sidx;
            if(isMinRep(repI,norm)){
                indexList.push_back(repI);
                locDim++;
                #ifdef KEEP_BASIS_NORM
                    normList.push_back(norm);
                #endif
            }
        }
    }
}

// construct subspace basis from reps loaded from file
void Basis::gen(std::string basisfile){
    subDim = 0;
    indexList.clear();
    read<ind_int>(&indexList, basisfile);
    subDim = indexList.size();
}

// construct subspace basis and norm from reps loaded from file
void Basis::gen(std::string basisfile, std::string normfile){
    subDim = 0;
    indexList.clear();
    read<ind_int>(&indexList, basisfile);
    subDim = indexList.size();
    
    normList.clear();
    read<double>(&normList, normfile);
    assert(normList.size()==indexList.size());
}
void Basis::gen(std::string basisfile, std::string normfile, int workerID, int workerNum){
    indexList.clear();
    read<ind_int>(&indexList, basisfile, workerID, workerNum);
    locDim = indexList.size();
    MPI_Allreduce(&locDim, &subDim, 1);
    repIStartList.resize(workerNum);
    repIEndList.resize(workerNum);
    MPI_Allgather(&indexList.front(), repIStartList.data(),1);
    MPI_Allgather(&indexList.back(), repIEndList.data(), 1);
    normList.clear();
    read<double>(&normList, normfile, workerID, workerNum);
    assert(normList.size()==indexList.size());
}

void Basis::saveBasis(std::string basisfile, bool is_app) {
    std::ofstream outfile;
    ind_int *d_pt = indexList.data();
    ind_int size = indexList.size();
    save<ind_int>(d_pt, size, &outfile, basisfile, is_app);
}

void Basis::saveBasis(std::string basisfile, std::string normfile, bool is_app) {
    std::ofstream outfile;
    save<ind_int>(indexList.data(), (ind_int)indexList.size(), &outfile, basisfile, is_app);
    save<double>(normList.data(), (ind_int)normList.size(), &outfile, normfile, is_app);
}

ind_int Basis::vecToRep(VecI& v) const {
    ind_int result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < N; i++){if (v.at(i) != 0) result += v[i] * mul[i];}
    return result;
}

pairIndex Basis::vecToRep(VecI& v, VecI& vp) const {
    ind_int r=0, rp=0;
    #pragma omp parallel for reduction(+:r,rp)
    for (int i = 0; i < N; i++){
        if (v.at(i) != 0) r += mul[i];
        if (vp.at(i) != 0) rp += mul[i];
    }
    pairIndex pairInd(r,rp);
    return pairInd;
}

void Basis::repToVec(ind_int index, VecI& v) const {
    for (int i = N - 1; i >= 0; i--){
        if (index > 0){
            v.at(i) = index % siteDim;
            index /= siteDim;
        }
        else{
            v.at(i) = 0;
        }
    }
}

void Basis::repToVec(pairIndex pairInd, VecI& v, VecI& vp) const {
    repToVec(pairInd.first,v);
    repToVec(pairInd.second,vp);
}

void Basis::repToVec(ind_int index, VecI& v, VecI& vp) const {
    repToVec(std::make_pair(fIndexList.at(index/sDim), sIndexList.at(index%sDim)), v, vp);
}

bool Basis::getBid(ind_int repI, int &bid) const {
    auto low = std::lower_bound(repIStartList.begin(),repIStartList.end(),repI);
    bid = low - repIStartList.begin();
    if(low!=repIStartList.end() && repI==repIStartList[bid]) return true;
    if(bid==0) return false;
    bid--;
    if(repI<=repIEndList[bid])return true;
    return false;
}

bool Basis::search(ind_int repI, ind_int &idx, const std::vector<ind_int> &indList) const {
    auto low = std::lower_bound(indList.begin(),indList.end(),repI);
    if((low==indList.end()) or (repI!=(*low))) return false;
    idx = low - indList.begin();
    return true;
}

bool Basis::search(ind_int repI, ind_int &idx) const {
    return search(repI, idx, indexList);
}

bool Basis::search(pairIndex pairRepI, ind_int &idx) const {
    assert(model != LATTICE_MODEL::HEISENBERG);
    auto fiter = fRepIdxHash.find(pairRepI.first); if(fiter==fRepIdxHash.end()) return false;
    auto siter = sRepIdxHash.find(pairRepI.second); if(fiter==sRepIdxHash.end()) return false;
    ind_int repI = (*fiter).second*sDim+(*siter).second;
    if (model==LATTICE_MODEL::HUBBARD and kIndex==-1){
        idx = repI;
        return true;
    }else{
        return search(repI, idx, indexList);
    }
    return false;
}

ind_int Basis::search(ind_int repI, const std::vector<ind_int> &indList) const {
    ind_int idx;
    if(search(repI,idx,indList)){
        return idx;
    } else{
        std::cerr<<"index not found in indList!"<<std::endl;
        exit(EXIT_FAILURE);
    }
}
ind_int Basis::search(ind_int repI) const {
    return search(repI, indexList);
}
ind_int Basis::search(pairIndex pairRepI) const {
    ind_int repI = fRepIdxHash.at(pairRepI.first)*sDim+sRepIdxHash.at(pairRepI.second);
    if(model==LATTICE_MODEL::HUBBARD){
        return repI;
    }else if(model==LATTICE_MODEL::t_J && !(pairRepI.first & pairRepI.second)){
        return repI;
    }else{
        std::cerr<<"pairRepI only defined for Hubbard and t_J!"<<std::endl;
        exit(EXIT_FAILURE);
    } 
}

/*
    *******************************************
    * Implementation for translation symmetry *
    *******************************************
*/
bool Basis::isMin(ind_int repI, VecI& symmList){
    symmList.clear();
    if(kIndex==-1) return true;
    if(PGRepIndex==-1){
        for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            ind_int repIp{0};
            for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
            }
            if (repI > repIp) return false;
            if (repI == repIp) symmList.push_back(r);
        }
    }else{
        for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                ind_int repIp{0};
                for (int i = 0; i < pt_lattice->getOrbNum();i++){
                    if(bitTest(repI,i)) bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                }
                if (repI > repIp) return false;
                if(repI == repI) {symmList.push_back(r); symmList.push_back(p);}
            }        
        }
    }
    return true;
}

bool Basis::isMinRep(ind_int repI, double& norm) const {
    // project out double occp
    if (model==LATTICE_MODEL::t_J){
        pairIndex pairRepI;
        pairRepI = getPairRepI(repI);
        if((pairRepI.first & pairRepI.second))return false;
    }
    // full hilbert space
    if (kIndex==-1) {
        norm = 1.0;
        return true;
    }
    // smallest index in the cycle?
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::t_J:{
            pairIndex pairRepI = getPairRepI(repI);
            if(PGRepIndex==-1){
                auto it = fMinRepSymmHash.find(pairRepI.first);
                if(it==fMinRepSymmHash.end()){
                    return false;
                }else{
                    for(auto symm=(*it).second.begin(); symm!=(*it).second.end(); symm++){
                        int r = *symm;
                        ind_int srepI=0;
                        for (int j = 0; j < pt_lattice->getOrbNum(); j++){
                            if(bitTest(pairRepI.second,j))bitSet(srepI,pt_lattice->getOrbTran(r,j));                            
                        }
                        if (pairRepI.second > srepI) return false;
                    }
                }
            }else{
                auto it = fMinRepSymmHash.find(pairRepI.first);
                if(it==fMinRepSymmHash.end()){
                    return false;
                }else{
                    for(auto symm=(*it).second.begin(); symm!=(*it).second.end(); symm += 2){
                        int r = *symm;
                        int p = *(symm+1);
                        ind_int srepI=0;
                        for (int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(pairRepI.second,i))bitSet(srepI,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));                      
                        }
                        std::cout<<"r:"<<r<<",p:"<<p<<", repI:"<<pairRepI.second<<", srepI:"<<srepI<<"\n";
                        if (pairRepI.second > srepI) return false;
                        std::cout<<"repI is in minRep:"<<repI<<", Hash size:"<<fMinRepSymmHash.size()<<"\n";
                    }
                }
            }

            // if(PGRepIndex==-1){
            //     for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            //         pairIndex pairRepIp{0,0};
            //         for (int i = 0; i < pt_lattice->getOrbNum(); i++){
            //             if(bitTest(pairRepI.first,i))bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));
            //         }
            //         if (pairRepI.first > pairRepIp.first) return false;
            //         else if (pairRepI.first == pairRepIp.first){
            //             for (int j = 0; j < pt_lattice->getOrbNum(); j++){
            //                 if(bitTest(pairRepI.second,j))bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,j));                            
            //             }
            //             if (pairRepI.second > pairRepIp.second) return false;
            //         }
            //     }
            // }else{
            //     for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            //         for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
            //             pairIndex pairRepIp{0,0};
            //             for (int i = 0; i < pt_lattice->getOrbNum();i++){
            //                 if(bitTest(pairRepI.first,i))bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
            //             }
            //             if (pairRepI.first > pairRepIp.first) return false;
            //             else if (pairRepI.first == pairRepIp.first){
            //                 for (int i = 0; i < pt_lattice->getOrbNum(); i++){
            //                     if(bitTest(pairRepI.second,i))bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));                      
            //                 }
            //                 if (pairRepI.second > pairRepIp.second) return false;
            //             }
            //         }    
            //     }
            // }
            break;
        }
        case LATTICE_MODEL::HEISENBERG:{
            #ifdef BINARY_REP
                if(PGRepIndex==-1){
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        ind_int repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                        }
                        if (repI > repIp) return false;
                    }
                }else{
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                            ind_int repIp{0};
                            for (int i = 0; i < pt_lattice->getOrbNum();i++){
                                if(bitTest(repI,i)) bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                            }
                            if (repI > repIp) return false;
                        }        
                    }
                }
            #else
                VecI initVec(N), finalVec(N), finalVec1(N);
                repToVec(repI, initVec);
                if(PGRepIndex==-1){
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        if (repI > vecToRep(finalVec)) return false;
                    }
                }else{
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                            for (int i = 0; i < pt_lattice->getOrbNum();i++){
                                finalVec1.at(pt_lattice->getOrbPG(p,i)) = finalVec.at(i);
                            }
                            if (repI > vecToRep(finalVec1)) return false;
                        }        
                    }
                }
            #endif

            break;
        }
    }
    // norm > 0?
    norm = (model==LATTICE_MODEL::HEISENBERG)?Norm(repI):minNorm(repI);
    std::cout<<"norm:"<<norm<<"\n";
    // sqrt(infinitismal)>>infinitesimal
    if (std::real(norm)*std::real(norm)>INFINITESIMAL){
        return true;
    }
    return false;
}

double Basis::minNorm(ind_int repI) const {
    if (kIndex==-1) return 1.0;
    if (model==LATTICE_MODEL::HUBBARD or model==LATTICE_MODEL::t_J){
        cdouble norm = 0.0;
        pairIndex pairRepI = getPairRepI(repI);
        auto repSymmIt = fMinRepSymmHash.find(pairRepI.first);
        assert_msg(repSymmIt!=fMinRepSymmHash.end(),"Basis::minNorm only defined for minimum repI in a cycle!");
        VecI seq, seqp;
        if(PGRepIndex==-1){
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm++){
                int r = *symm;
                pairIndex pairRepIp{0,0};
                seq.clear();
                seqp.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                    if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                    if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));}
                }
                if (pairRepIp==pairRepI) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r);
            }
            norm /= pt_lattice->getSiteNum();
        }else{
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm += 2){
                int r = *symm, p = *(symm+1);
                pairIndex pairRepIp{0,0};
                seq.clear();
                seqp.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                    if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                    if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                } 
                if (pairRepIp==pairRepI) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
            }
            norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
        }
        assert(std::abs(std::imag(norm))<INFINITESIMAL);
        return std::sqrt(std::real(norm));
    }else{
        std::cout<<"Basis::minNorm only defined fro Hubbard and t_J"<<std::endl;
        exit(1);
    }
}
double Basis::Norm(ind_int repI) const {
    if (kIndex==-1) return 1.0;
    cdouble norm = 0.0;
    switch(model){
        case HUBBARD:case t_J:{
            pairIndex pairRepI = getPairRepI(repI);
            VecI seq, seqp;
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIndex pairRepIp{0,0};
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                        if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));}
                    }
                    if (pairRepIp==pairRepI) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIndex pairRepIp{0,0};
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        if (pairRepIp==pairRepI) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }
                }
                norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
            break;
        }

        case HEISENBERG:{
        #ifdef BINARY_REP
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    ind_int repIp{0};
                    for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                        if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                    }
                    if (repIp==repI) norm += pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        ind_int repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                        }
                        if (repIp==repI) norm += pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }
                }
                norm /= pt_lattice->getSiteNum()* pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
        #else
            VecI initVec(N), finalVec(N);
            repToVec(repI, initVec);
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int i = 0; i < N; i++){
                        finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                    }
                    if (vecToRep(finalVec)==repI) norm += pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        for(int i = 0; i < N; i++){
                            finalVec.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVec.at(i);
                        }
                        if (vecToRep(finalVec)==repI) norm += pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }   
                }
                norm /= pt_lattice->getSiteNum()* pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
        #endif
            break;
        }
    }
    std::cout<<"Undefined Model!"<<std::endl;
    exit(1);
}

// For a given ind, apply translation symmetry, return all the resulting indexes in finalInd
void Basis::genSymm(ind_int ind, std::vector<ind_int>& finalInd) const {
    /*
        ind is repInt for Heisenberg model
        ind is the idx for pairRepI of Hubbard model
    */
    if (kIndex == -1){finalInd.push_back(ind);return;}
    switch(model){
        case HUBBARD:case t_J:{
            pairIndex repI;
            repI = getPairRepI(ind);
            for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                pairIndex repIp{0,0};
                for(int i = 0; i < N; i++){
                    if(bitTest(repI.first,i)) bitSet(repIp.first,pt_lattice->getOrbTran(r,i));
                    if(bitTest(repI.second,i)) bitSet(repIp.second,pt_lattice->getOrbTran(r,i));
                }
                finalInd.push_back(search(repIp));
            }
            break;
        }

        case HEISENBERG:{
            #ifdef BINARY_REP
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    ind_int indp{0};
                    for(int i = 0; i<N; i++){
                        if(bitTest(ind,i))bitSet(indp,pt_lattice->getOrbTran(r,i));
                    }
                    finalInd.push_back(indp);
                }
            
            #else
                VecI veci(N), vecf(N);
                repToVec(ind, veci);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int i = 0; i<N; i++){
                        vecf.at(pt_lattice->getOrbTran(r,i)) = veci.at(i);
                    }
                    finalInd.push_back(vecToRep(vecf));
                }
            #endif

            break;
        }

        default:{
            std::cout<<"model not defined! must be: HUBBARD,t_J,HEISENBERG."<<std::endl;
            exit(1);
            break;
        }
    }
}

void Basis::genSymm(ind_int rowid, std::vector<ind_int>& finalInd, std::vector<cdouble>& factorList) const {
    ind_int repI = getRepI(rowid);
    if (kIndex == -1){finalInd.push_back(repI); factorList.push_back(1.0);return;}
    cdouble initNorm = getNorm(rowid);
    switch(model){
        case LATTICE_MODEL::HUBBARD:{
            pairIndex pairRepI;
            pairRepI = getPairRepI(repI);
            VecI seq, seqp;
            if(PGRepIndex==-1){
                initNorm *= pt_lattice->getSiteNum();
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIndex pairRepIp{0,0};
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                        if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));}
                    }
                    finalInd.push_back(search(pairRepIp));
                    factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm*seqSign(seq) * seqSign(seqp));
                }
            }else{
                initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIndex pairRepIp{0,0};
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        finalInd.push_back(search(pairRepIp));
                        factorList.push_back(pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p)/initNorm*seqSign(seq) * seqSign(seqp));
                    }
                }
            }
            break;
        }

        case LATTICE_MODEL::t_J:{
            pairIndex pairRepI;
            pairRepI = getPairRepI(repI);
            VecI seq;
            if(PGRepIndex==-1){
                initNorm *= pt_lattice->getSiteNum();
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIndex pairRepIp{0,0};
                    seq.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                        else if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                    }
                    finalInd.push_back(search(pairRepIp));
                    factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm*seqSign(seq));
                }
            }else{
                initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIndex pairRepIp{0,0};
                        seq.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            else if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        finalInd.push_back(search(pairRepIp));
                        factorList.push_back(pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p)/initNorm*seqSign(seq));
                    }
                }
            }
            break;
        }

        case LATTICE_MODEL::HEISENBERG:{
            #ifdef BINARY_REP
                if(PGRepIndex==-1){
                    initNorm *= pt_lattice->getSiteNum();
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        ind_int repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                        }
                        finalInd.push_back(repIp);
                        factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm);
                    }
                }else{
                    initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                            ind_int repIp{0};
                            for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                                if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                            }
                            finalInd.push_back(repIp);
                            factorList.push_back(pt_lattice->expKR(kIndex,r)* pt_lattice->getChi(PGRepIndex,p)/initNorm);
                        }
                    }
                }
            
            #else
                VecI initVec(N), finalVec(N);
                repToVec(repI, initVec);
                if(PGRepIndex==-1){
                    initNorm *= pt_lattice->getSiteNum();
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < N; i++){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        finalInd.push_back(vecToRep(finalVec));
                        factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm);
                    }
                }else{
                    initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                            for(int i = 0; i < N; i++){
                                finalVec.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVec.at(i);
                            }
                            finalInd.push_back(vecToRep(finalVec));
                            factorList.push_back(pt_lattice->expKR(kIndex,r)* pt_lattice->getChi(PGRepIndex,p)/initNorm);
                        }
                    }
                }
            #endif

            break;
        }
        
        default:{
            std::cout<<"model not defined! must be: HUBBARD,t_J,HEISENBERG."<<std::endl;
            exit(1);
            break;
        }
    }
}