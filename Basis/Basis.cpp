//
//  Basis.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Basis.hpp"

// distribute workload among mpi workers
void work_load(idx_t size, int workerID, int workerNum, idx_t& idxStart, idx_t& idxEnd){
    assert((workerNum>0) && (workerID<workerNum));
    idx_t nlocmax = (size + workerNum - 1)/workerNum;
    idxStart = workerID * nlocmax;
    idxEnd = (idxStart + nlocmax)<size?(idxStart + nlocmax):size;
}

/*
    ***************
    * Basis Class *
    ***************
*/
// initialize
Basis::Basis(LATTICE_MODEL input_model, Geometry *pt_lat, VecI occList, int kInd, int PGRepInd, int siteD):\
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
        case LATTICE_MODEL::HUBBARD: {
            Nocc = occList;
            assert(Nocc.at(0)<=N);
            assert(Nocc.at(1)<=N);
            fDim = combination<idx_t>((idx_t)N, (idx_t)(Nocc[0]));
            sDim = combination<idx_t>((idx_t)N, (idx_t)(Nocc[1]));
            totDim = fDim * sDim;
            break;
        }
        case LATTICE_MODEL::tJ: {
            Nocc = occList;
            assert((Nocc.at(0)+Nocc.at(1))<=N);
            fDim = combination<idx_t>((idx_t)N, (idx_t)(Nocc[0]));
            sDim = combination<idx_t>((idx_t)N, (idx_t)(Nocc[1]));
            totDim = fDim * combination<idx_t>((idx_t)(N-Nocc[0]), (idx_t)(Nocc[1]));
            break;
        }
        case LATTICE_MODEL::HEISENBERG: {
            Sztot = occList.at(0) - occList.at(1);
            Nocc = occList;
            #ifdef BINARY_REP
                assert((Nocc.at(0) + Nocc.at(1))==N);
                // calculate totDim
                totDim = 1;
                idx_t n = N;
                for (int i = 0; i < siteDim; ++i){
                    totDim *= combination<idx_t>(n, (idx_t)(Nocc[i]));
                    n -= Nocc[i];
                }
            #else
                // initialize vec
                for (idx_t i = 0; i < siteDim; ++i){
                    for (int j = 0; j < occList[i]; j++) vec.push_back(i);
                }
            #endif
            break;
        }
        default:
            break;
    }
    initMinMaxRep();
    // default full hilbert space
    subDim = totDim;
    locDim = subDim;
    // initialize mul
    for (int i = 0; i < N; ++i) mul.push_back(pow((idx_t) siteDim, (idx_t) (N - i - 1)));  
    // generate fIndexList and sIndexList for Hubbard and tJ
    gendcmp(); 
}
void Basis::initMinMaxRep() const {
#ifdef BINARY_REP
    assert(siteDim==2);
    if(model==LATTICE_MODEL::HUBBARD or model==LATTICE_MODEL::tJ){
        fminRep=0;fmaxRep=0;
        for(int i=0; i<Nocc.at(0); ++i){
            bitSet(fminRep,i);
            bitSet(fmaxRep,N-1-i);
        }
        sminRep=0;smaxRep=0;
        for(int i=0; i<Nocc.at(1); ++i){
            bitSet(sminRep,i);
            bitSet(smaxRep,N-1-i);
        }
    }
    if(model==LATTICE_MODEL::HEISENBERG){
        fminRep=0;fmaxRep=0;
        for(int i=0; i<Nocc.at(1); ++i){
            bitSet(fminRep,i);
            bitSet(fmaxRep,N-1-i);
        }
    }
#else
    vec.clear();
    vecp.clear();
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::tJ:
            for (idx_t i = 0; i < (N - Nocc); ++i) vec.push_back(0);
            for (idx_t i = 0; i < Nocc; ++i) vec.push_back(1);
            for (idx_t i = 0; i < (Np - Npocc); ++i) vecp.push_back(0);
            for (idx_t i = 0; i < Npocc; ++i) vecp.push_back(1);
            break;
        default:break;
    }
#endif
}

void Basis::gendcmp(){
    if(model==LATTICE_MODEL::HUBBARD || model==LATTICE_MODEL::tJ){
        fIndexList.clear(); sIndexList.clear();
        fIndexList.reserve(fDim); sIndexList.reserve(sDim);

        idx_t counter = 0;
        idx_t repI = fminRep;
        while(repI<=fmaxRep){
            fIndexList.push_back(repI);
            fRepIdxHash[repI] = counter;

            repI = nextLexicographicalNumber<idx_t>(repI);
            counter++;
        }
        assert(counter == fDim);

        counter = 0;
        repI = sminRep;
        while(repI<=smaxRep){
            sIndexList.push_back(repI);
            sRepIdxHash[repI] = counter;

            repI = nextLexicographicalNumber<idx_t>(repI);
            counter++;
        }
        assert(counter == sDim);

        for(idx_t idx=0; idx<fDim; idx++){
            idx_t repI = fIndexList.at(idx);
            VecI symmList;
            if(isMin(repI, symmList)){
                fMinRepSymmHash[repI] = symmList;
                fMinRepList.push_back(repI);
            } 
        }
    }
}
// generate Basis for the subspace labeled by kInd
void Basis::construct(){
    initMinMaxRep();
    subDim = 0;
    indexList.clear();
    normList.clear();
    idx_t repI;
    double norm;
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::tJ:{
            if (!(kIndex==-1 and model==LATTICE_MODEL::HUBBARD)){ 
                for(idx_t fidx=0;fidx<fDim;fidx++){
                    if(!isfMin(fIndexList[fidx])) continue;
                    for(idx_t sidx=0;sidx<sDim;sidx++){
                        repI = fidx*sDim+sidx;
                        if(isMinRep(repI,norm)){
                            // std::cout<<"test repI is minRepI"<<"\n";
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
        case LATTICE_MODEL::HEISENBERG:{
        #ifdef BINARY_REP
            idx_t repI = fminRep;
            while(repI<=fmaxRep){
                if (isMinRep(repI, norm)){
                    indexList.push_back(repI);
                    subDim++;
                    #ifdef KEEP_BASIS_NORM
                        if (kIndex!=-1) normList.push_back(norm);
                    #endif
                }
                repI = nextLexicographicalNumber<idx_t>(repI);
            }
        #else
            // do{
            //     idx_t repI = vecToRep(vec);
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

void Basis::construct(int workerID, int workerNum){
    locDim = 0;
    assert(model==LATTICE_MODEL::HUBBARD or model==LATTICE_MODEL::tJ);
    if(model==LATTICE_MODEL::HUBBARD) assert(kIndex!=-1);
    idx_t repI;
    double norm;
    idx_t size = fMinRepList.size();
    idx_t fidxStart, fidxEnd;
    work_load(size, workerID, workerNum, fidxStart, fidxEnd);
    for(idx_t idx=fidxStart; idx<fidxEnd; idx++){
        idx_t fidx = fRepIdxHash.at(fMinRepList[idx]);
        for(idx_t sidx=0;sidx<sDim;sidx++){
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
void Basis::construct(std::string basisfile){
    subDim = 0;
    indexList.clear();
    read<idx_t>(&indexList, basisfile);
    subDim = indexList.size();
}

// construct subspace basis and norm from reps loaded from file
void Basis::construct(std::string basisfile, std::string normfile){
    subDim = 0;
    indexList.clear();
    read<idx_t>(&indexList, basisfile);
    subDim = indexList.size();
    
    normList.clear();
    read<double>(&normList, normfile);
    assert(normList.size()==indexList.size());
}
void Basis::construct(std::string basisfile, std::string normfile, int workerID, int workerNum){
    indexList.clear();
    read<idx_t>(&indexList, basisfile, workerID, workerNum);
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

void Basis::construct(bool saved, std::string basisDir) {
    if (saved) {
        construct(basisDir + "/basis", basisDir + "/norm");
    } else {
        construct();
    }
}



void Basis::print(std::ostream& os) const {
    os<<"Basis Info:\n"<<model<<"\nkid:"<<kIndex<<". pid:"<<PGRepIndex<<"\n";
    os<<"tot/sub/loc dim:"<<totDim<<" / "<<subDim<<" / "<<locDim<<"\n";
}

void Basis::save(std::string basisDir, bool saveNorm, bool is_app) {
    mkdir_fs(basisDir);
    std::ofstream outfile;
    ::save<idx_t>(indexList.data(), (idx_t)indexList.size(), &outfile, basisDir + "/basis", is_app);
    if (saveNorm) {
        ::save<double>(normList.data(), (idx_t)normList.size(), &outfile, basisDir + "/norm", is_app); 
    }
}

idx_t Basis::vecToRep(VecI& v) const {
    idx_t result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < N; ++i){if (v.at(i) != 0) result += v[i] * mul[i];}
    return result;
}

pairIdx_t Basis::vecToRep(VecI& v, VecI& vp) const {
    idx_t r=0, rp=0;
    #pragma omp parallel for reduction(+:r,rp)
    for (int i = 0; i < N; ++i){
        if (v.at(i) != 0) r += mul[i];
        if (vp.at(i) != 0) rp += mul[i];
    }
    pairIdx_t pairInd(r,rp);
    return pairInd;
}

void Basis::repToVec(idx_t index, VecI& v) const {
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

void Basis::repToVec(pairIdx_t pairInd, VecI& v, VecI& vp) const {
    repToVec(pairInd.first,v);
    repToVec(pairInd.second,vp);
}

void Basis::repToVec(idx_t index, VecI& v, VecI& vp) const {
    repToVec(std::make_pair(fIndexList.at(index/sDim), sIndexList.at(index%sDim)), v, vp);
}

bool Basis::getBid(idx_t repI, int &bid) const {
    auto low = std::lower_bound(repIStartList.begin(),repIStartList.end(),repI);
    bid = low - repIStartList.begin();
    if(low!=repIStartList.end() && repI==repIStartList[bid]) return true;
    if(bid==0) return false;
    bid--;
    if(repI<=repIEndList[bid])return true;
    return false;
}

bool Basis::search(idx_t repI, idx_t &idx, const std::vector<idx_t> &indList) const {
    auto low = std::lower_bound(indList.begin(),indList.end(),repI);
    if((low==indList.end()) or (repI!=(*low))) return false;
    idx = low - indList.begin();
    return true;
}

bool Basis::search(idx_t repI, idx_t &idx) const {
    return search(repI, idx, indexList);
}

bool Basis::search(pairIdx_t pairRepI, idx_t &idx) const {
    assert(model != LATTICE_MODEL::HEISENBERG);
    auto fiter = fRepIdxHash.find(pairRepI.first); if(fiter==fRepIdxHash.end()) return false;
    auto siter = sRepIdxHash.find(pairRepI.second); if(fiter==sRepIdxHash.end()) return false;
    idx_t repI = (*fiter).second*sDim+(*siter).second;
    if (model==LATTICE_MODEL::HUBBARD and kIndex==-1){
        idx = repI;
        return true;
    }else{
        return search(repI, idx, indexList);
    }
    return false;
}

idx_t Basis::search(idx_t repI, const std::vector<idx_t> &indList) const {
    idx_t idx;
    if(search(repI,idx,indList)){
        return idx;
    } else{
        std::cerr<<"index not found in indList!"<<std::endl;
        exit(EXIT_FAILURE);
    }
}
idx_t Basis::search(idx_t repI) const {
    return search(repI, indexList);
}
idx_t Basis::search(pairIdx_t pairRepI) const {
    idx_t repI = fRepIdxHash.at(pairRepI.first)*sDim+sRepIdxHash.at(pairRepI.second);
    if(model==LATTICE_MODEL::HUBBARD){
        return repI;
    }else if(model==LATTICE_MODEL::tJ && !(pairRepI.first & pairRepI.second)){
        return repI;
    }else{
        std::cerr<<"pairRepI only defined for Hubbard and tJ!"<<std::endl;
        exit(EXIT_FAILURE);
    } 
}

/*
    *******************************************
    * Implementation for translation symmetry *
    *******************************************
*/
bool Basis::isMin(idx_t repI, VecI& symmList){
    symmList.clear();
    if(kIndex==-1) return true;
    if(PGRepIndex==-1){
        for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            idx_t repIp{0};
            for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
            }
            if (repI > repIp) return false;
            if (repI == repIp) symmList.push_back(r);
        }
    }else{
        for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                idx_t repIp{0};
                for (int i = 0; i < pt_lattice->getOrbNum();++i){
                    if(bitTest(repI,i)) bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                }
                if (repI > repIp) return false;
                if(repI == repIp) {symmList.push_back(r); symmList.push_back(p);}
            }        
        }
    }
    return true;
}

bool Basis::isMinRep(idx_t repI, double& norm) const {
    // project out double occp
    if (model==LATTICE_MODEL::tJ){
        pairIdx_t pairRepI;
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
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI = getPairRepI(repI);
            if(PGRepIndex==-1){
                auto it = fMinRepSymmHash.find(pairRepI.first);
                if(it==fMinRepSymmHash.end()){
                    return false;
                }else{
                    for(auto symm=(*it).second.begin(); symm!=(*it).second.end(); symm++){
                        int r = *symm;
                        idx_t srepI=0;
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
                        idx_t srepI=0;
                        for (int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(pairRepI.second,i))bitSet(srepI,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));                      
                        }
                        if (pairRepI.second > srepI) return false;
                    }
                }
            }

            // if(PGRepIndex==-1){
            //     for (int r = 0; r < pt_lattice->getSiteNum(); r++){
            //         pairIdx_t pairRepIp{0,0};
            //         for (int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
            //             pairIdx_t pairRepIp{0,0};
            //             for (int i = 0; i < pt_lattice->getOrbNum();++i){
            //                 if(bitTest(pairRepI.first,i))bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
            //             }
            //             if (pairRepI.first > pairRepIp.first) return false;
            //             else if (pairRepI.first == pairRepIp.first){
            //                 for (int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
                        idx_t repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                        }
                        if (repI > repIp) return false;
                    }
                }else{
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                            idx_t repIp{0};
                            for (int i = 0; i < pt_lattice->getOrbNum();++i){
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
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        if (repI > vecToRep(finalVec)) return false;
                    }
                }else{
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                            for (int i = 0; i < pt_lattice->getOrbNum();++i){
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
    // std::cout<<"norm:"<<norm<<"\n";
    // sqrt(infinitismal)>>infinitesimal
    if (std::real(norm)*std::real(norm)>INFINITESIMAL){
        return true;
    }
    return false;
}

double Basis::minNorm(idx_t repI) const {
    if (kIndex==-1) return 1.0;
    cdouble norm = 0.0;
    pairIdx_t pairRepI = getPairRepI(repI);
    auto repSymmIt = fMinRepSymmHash.find(pairRepI.first);
    assert_msg(repSymmIt!=fMinRepSymmHash.end(),"Basis::minNorm only defined for minimum repI in a cycle!");
    if (model==LATTICE_MODEL::HUBBARD){
        VecI seq, seqp;
        if(PGRepIndex==-1){
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm++){
                int r = *symm;
                pairIdx_t pairRepIp{0,0};
                cdouble tbc=1.0;
                seq.clear();
                seqp.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                    if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                    if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                }
                if (pairRepIp==pairRepI) norm += tbc * seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r);
            }
            norm /= pt_lattice->getSiteNum();
        }else{
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm += 2){
                int r = *symm, p = *(symm+1);
                pairIdx_t pairRepIp{0,0};
                seq.clear();
                seqp.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                    if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                    if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                } 
                if (pairRepIp==pairRepI) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
            }
            norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
        }
        assert(std::abs(std::imag(norm))<INFINITESIMAL);
        return std::sqrt(std::real(norm));
    }else if(model==LATTICE_MODEL::tJ){
        VecI seq;
        if(PGRepIndex==-1){
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm++){
                int r = *symm;
                pairIdx_t pairRepIp{0,0};
                seq.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                    if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                    else if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                }
                if (pairRepIp==pairRepI) norm += seqSign(seq) * pt_lattice->expKR(kIndex,r);
            }
            norm /= pt_lattice->getSiteNum();
        }else{
            for (auto symm = (*repSymmIt).second.begin(); symm != (*repSymmIt).second.end(); symm += 2){
                int r = *symm, p = *(symm+1);
                pairIdx_t pairRepIp{0,0};
                seq.clear();
                for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                    if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                    else if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                } 
                if (pairRepIp==pairRepI) norm += seqSign(seq) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
            }
            norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
        }
        assert(std::abs(std::imag(norm))<INFINITESIMAL);
        return std::sqrt(std::real(norm));
    }else{
        std::cout<<"Basis::minNorm only defined fro Hubbard and tJ"<<std::endl;
        exit(1);
    }
}
double Basis::Norm(idx_t repI) const {
    if (kIndex==-1) return 1.0;
    cdouble norm = 0.0;
    switch(model){
        case LATTICE_MODEL::HUBBARD:{
            pairIdx_t pairRepI = getPairRepI(repI);
            VecI seq, seqp;
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    cdouble tbc=1.0;
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                        if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                    }
                    if (pairRepIp==pairRepI) norm += tbc * seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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

        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI = getPairRepI(repI);
            VecI seq;
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    seq.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                        else if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                    }
                    if (pairRepIp==pairRepI) norm += seqSign(seq) * pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            else if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        if (pairRepIp==pairRepI) norm += seqSign(seq) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }
                }
                norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
            break;
        }

        case LATTICE_MODEL::HEISENBERG:{
        #ifdef BINARY_REP
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    idx_t repIp{0};
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                    }
                    if (repIp==repI) norm += pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        idx_t repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
                    for(int i = 0; i < N; ++i){
                        finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                    }
                    if (vecToRep(finalVec)==repI) norm += pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        for(int i = 0; i < N; ++i){
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
void Basis::genSymm(idx_t ind, std::vector<idx_t>& finalInd) const {
    /*
        ind is repInt for Heisenberg model
        ind is the idx for pairRepI of Hubbard model
    */
    if (kIndex == -1){finalInd.push_back(ind);return;}
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::tJ:{
            pairIdx_t repI;
            repI = getPairRepI(ind);
            for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                pairIdx_t repIp{0,0};
                for(int i = 0; i < N; ++i){
                    if(bitTest(repI.first,i)) bitSet(repIp.first,pt_lattice->getOrbTran(r,i));
                    if(bitTest(repI.second,i)) bitSet(repIp.second,pt_lattice->getOrbTran(r,i));
                }
                finalInd.push_back(search(repIp));
            }
            break;
        }

        case LATTICE_MODEL::HEISENBERG:{
            #ifdef BINARY_REP
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    idx_t indp{0};
                    for(int i = 0; i<N; ++i){
                        if(bitTest(ind,i))bitSet(indp,pt_lattice->getOrbTran(r,i));
                    }
                    finalInd.push_back(indp);
                }
            
            #else
                VecI veci(N), vecf(N);
                repToVec(ind, veci);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int i = 0; i<N; ++i){
                        vecf.at(pt_lattice->getOrbTran(r,i)) = veci.at(i);
                    }
                    finalInd.push_back(vecToRep(vecf));
                }
            #endif

            break;
        }

        default:{
            std::cout<<"model not defined! must be: HUBBARD,tJ,HEISENBERG."<<std::endl;
            exit(1);
            break;
        }
    }
}

void Basis::genSymm(idx_t rowid, std::vector<idx_t>& finalInd, std::vector<cdouble>& factorList) const {
    idx_t repI = getRepI(rowid);
    if (kIndex == -1){finalInd.push_back(repI); factorList.push_back(1.0);return;}
    cdouble initNorm = getNorm(rowid);
    switch(model){
        case LATTICE_MODEL::HUBBARD:{
            pairIdx_t pairRepI;
            pairRepI = getPairRepI(repI);
            VecI seq, seqp;
            if(PGRepIndex==-1){
                initNorm *= pt_lattice->getSiteNum();
                cdouble tbc = 1.0;
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                        if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                    }
                    finalInd.push_back(search(pairRepIp));
                    factorList.push_back(tbc*pt_lattice->expKR(kIndex,r)/initNorm*seqSign(seq) * seqSign(seqp));
                }
            }else{
                initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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

        case LATTICE_MODEL::tJ:{
            pairIdx_t pairRepI;
            pairRepI = getPairRepI(repI);
            VecI seq;
            if(PGRepIndex==-1){
                initNorm *= pt_lattice->getSiteNum();
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    seq.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
                        idx_t repIp{0};
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(repI,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                        }
                        finalInd.push_back(repIp);
                        factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm);
                    }
                }else{
                    initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                            idx_t repIp{0};
                            for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
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
                        for(int i = 0; i < N; ++i){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        finalInd.push_back(vecToRep(finalVec));
                        factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm);
                    }
                }else{
                    initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                            for(int i = 0; i < N; ++i){
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
            std::cout<<"model not defined! must be: HUBBARD,tJ,HEISENBERG."<<std::endl;
            exit(1);
            break;
        }
    }
}

void Basis::genSymm(idx_t rowid, std::vector<pairIdx_t>& finalInd, std::vector<cdouble>& factorList) const {
    idx_t repI = getRepI(rowid);
    auto pairRepI = getPairRepI(repI);
    if (kIndex == -1) {
        finalInd.push_back(pairRepI); 
        factorList.push_back(1.0);
        return;
    }
    cdouble initNorm = getNorm(rowid);
    switch (model) {
        case LATTICE_MODEL::HUBBARD: {
            VecI seq, seqp;
            if (PGRepIndex==-1) {
                initNorm *= pt_lattice->getSiteNum();
                cdouble tbc = 1.0;
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                        if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seqp.push_back(pt_lattice->getOrbTran(r,i));tbc *= pt_lattice->getOrbTranPhase(r,i);}
                    }
                    finalInd.push_back(pairRepIp);
                    factorList.push_back(tbc*pt_lattice->expKR(kIndex,r)/initNorm*seqSign(seq) * seqSign(seqp));
                }
            }else{
                initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        finalInd.push_back(pairRepIp);
                        factorList.push_back(pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p)/initNorm*seqSign(seq) * seqSign(seqp));
                    }
                }
            }
            break;
        }

        case LATTICE_MODEL::tJ:{
            VecI seq;
            if(PGRepIndex==-1){
                initNorm *= pt_lattice->getSiteNum();
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIdx_t pairRepIp{0,0};
                    seq.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                        if(bitTest(pairRepI.first,i)) {bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                        else if(bitTest(pairRepI.second,i)) {bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,i));seq.push_back(pt_lattice->getOrbTran(r,i));}
                    }
                    finalInd.push_back(pairRepIp);
                    factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm*seqSign(seq));
                }
            }else{
                initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        pairIdx_t pairRepIp{0,0};
                        seq.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); ++i){
                            if(bitTest(pairRepI.first,i)){bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                            else if(bitTest(pairRepI.second,i)){bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));}
                        } 
                        finalInd.push_back(pairRepIp);
                        factorList.push_back(pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p)/initNorm*seqSign(seq));
                    }
                }
            }
            break;
        }
        
        default:{
            std::cout<<"pairIdx_t not defiend for "<<model<<"\n";
            exit(1);
            break;
        }
    }
}