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
    N = pt_lattice->getOrbNum();
    assert(kIndex<pt_lattice->getSiteNum());
    assert(PGRepIndex<pt_lattice->getPGRepNum());
    Np = N;
    Nocc = occList.at(0);
    Npocc = occList.at(1);
    switch(model){
        case LATTICE_MODEL::HUBBARD:
            assert(Nocc<=N);
            assert(Npocc<=Np);
            fDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc));
            sDim = combination<ind_int>((ind_int)Np, (ind_int)(Npocc));
            totDim = fDim * sDim;
            break;
        case LATTICE_MODEL::t_J:
            assert((Nocc+Npocc)<=N);
            fDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc));
            sDim = combination<ind_int>((ind_int)N, (ind_int)(occList[1]));
            totDim = fDim * combination<ind_int>((ind_int)(N-Nocc), (ind_int)(occList[1]));
            break;
        // case t_J:
        //     assert((Nocc+occList[1]+occList[2])==N);
        //     Np = occList[1] + occList[2];
        //     for (ind_int i = 0; i < (N - Nocc); i++) vec.push_back(0);
        //     for (ind_int i = 0; i < Nocc; i++) vec.push_back(1);
        //     for (ind_int i = 0; i < (occList[1]); i++) vecp.push_back(0);
        //     for (ind_int i = 0; i < occList[2]; i++) vecp.push_back(1);
        //     fDim = combination<ind_int>((ind_int)N, (ind_int)(Nocc));
        //     sDim = combination<ind_int>((ind_int)(N-Nocc), (ind_int)(occList[1]));
        //     totDim = fDim * sDim;
        //     initVec.resize(N,0); initVecp.resize(N,0);
        //     finalVec.resize(N,0); finalVecp.resize(N,0);
        //     break;
        case LATTICE_MODEL::HEISENBERG:
            assert((Nocc+Npocc)==N);
            // initialize vec
            for (ind_int i = 0; i < siteDim; i++){
                for (int j = 0; j < occList[i]; j++) vec.push_back(i);
            }
            // calculate totDim
            totDim = 1;
            ind_int n = N;
            for (int i = 0; i < siteDim; i++){
                totDim *= combination<ind_int>(n, (ind_int)(occList[i]));
                n -= occList[i];
            }
            break;
        default:break;
    }
    initStartVec();
    gendcmp();
    // default full hilbert space
    subDim = totDim;
    locDim = subDim;
    // initialize mul
    for (int i = 0; i < N; i++) mul.push_back(pow((ind_int) siteDim, (ind_int) (N - i - 1)));   
}
void Basis::initStartVec() const {
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
}
// generate Basis for the full Hilbert Space
void Basis::genFull(){
    ind_int counter = 0;
    indexList.clear();
    fIndexList.clear();
    sIndexList.clear();
    switch(model){
        case HUBBARD:{
            fIndexList.reserve(fDim);
            sIndexList.reserve(sDim);
            counter = 0;
            do{
                fIndexList.push_back(vecToInd(vec));
                counter++;
            }while (std::next_permutation(vec.begin(), vec.end()));
            assert(counter == fDim);
            
            counter = 0;
            do{
                sIndexList.push_back(vecToInd(vecp));
                counter++;
            }while (std::next_permutation(vecp.begin(), vecp.end()));
            assert(counter == sDim);
            break;
        }
        
        case t_J:{
            fIndexList.reserve(fDim);
            indexList.reserve(totDim);
            sIndexList.reserve(sDim);
            counter = 0;
            do{
                fIndexList.push_back(vecToInd(vec));
                counter++;
            }while (std::next_permutation(vec.begin(), vec.end()));
            assert(counter == fDim);

            counter = 0;
            do{
                sIndexList.push_back(vecToInd(vecp));
                counter++;
            }while (std::next_permutation(vecp.begin(), vecp.end()));
            assert(counter == sDim);

            // project out double occupancy
            counter = 0;
            VecI initVec(N), initVecp(N);
            for (ind_int i = 0; i < fDim; i++){
                for (ind_int j = 0; j < sDim; j++){
                    if(fIndexList.at(i) & sIndexList.at(j)) continue;
                    indexList.push_back(i * sDim + j);
                    counter++;
                }
            }
            assert(totDim==counter);
            break;
        }
        case HEISENBERG:{
            indexList.reserve(totDim);
            counter = 0;
            do{
                indexList.push_back(vecToInd(vec));
                counter++;
            }while (std::next_permutation(vec.begin(), vec.end()));
            assert(counter == totDim);
            break;
        }
    }  
}
void Basis::gendcmp(){
    if(model==HUBBARD or model==t_J){
        fIndexList.clear();
        sIndexList.clear();
        fIndexList.reserve(fDim);
        sIndexList.reserve(sDim);
        ind_int counter = 0;
        do{
            fIndexList.push_back(vecToInd(vec));
            counter++;
        }while (std::next_permutation(vec.begin(), vec.end()));
        assert(counter == fDim);
        counter = 0;
        do{
            sIndexList.push_back(vecToInd(vecp));
            counter++;
        }while (std::next_permutation(vecp.begin(), vecp.end()));
        assert(counter == sDim);
    }
}
// generate Basis for the subspace labeled by kInd
void Basis::gen(){
    subDim = 0;
    indexList.clear();
    fIndexList.clear();
    sIndexList.clear();
    double norm;
    normList.clear();
    switch(model){
        case HUBBARD:case t_J:{
            fIndexList.reserve(fDim);
            sIndexList.reserve(sDim);
            ind_int counter = 0;
            do{
                fIndexList.push_back(vecToInd(vec));
                counter++;
            }while (std::next_permutation(vec.begin(), vec.end()));
            assert(counter == fDim);
            counter = 0;
            do{
                sIndexList.push_back(vecToInd(vecp));
                counter++;
            }while (std::next_permutation(vecp.begin(), vecp.end()));
            assert(counter == sDim);
            if (kIndex!=-1){
                initStartVec();
                pairIndex pairInd;
                int count = 0;
                do{
                    do{
                        pairInd = vecToInd(vec,vecp);
                        ind_int initInd = search(pairInd);
                        if (isRep(initInd, norm)){
                            indexList.push_back(initInd);
                            subDim++;
                            #ifdef KEEP_BASIS_NORM
                                normList.push_back(norm);
                            #endif
                        }
                    }while(std::next_permutation(vecp.begin(), vecp.end()));
                }while(std::next_permutation(vec.begin(), vec.end()));
            }else if (model==LATTICE_MODEL::t_J){
                // project out double occupancy
                counter = 0;
                VecI initVec(N), initVecp(N);
                for (ind_int i = 0; i < fDim; i++){
                    indToVec(fIndexList.at(i), initVec);
                    for (ind_int j = 0; j < sDim; j++){
                        indToVec(sIndexList.at(j), initVec);
                        bool no_double_occp = true;
                        for (int n = 0; n < N; n++){
                            if ((initVec.at(n) + initVecp.at(n)) > 1){
                                no_double_occp = false;
                                break;
                            }
                        }
                        if (no_double_occp){
                            indexList.push_back(i * sDim + j);
                            counter++;
                        } 
                    }
                }
                assert(totDim==counter);
            }
            break;
        }
        case HEISENBERG:{
            do{
                ind_int initInd = vecToInd(vec);
                if (isRep(initInd, norm)){
                    indexList.push_back(initInd);
                    subDim++;
                    #ifdef KEEP_BASIS_NORM
                        if (kIndex!=-1) normList.push_back(norm);
                    #endif
                }
            }while (std::next_permutation(vec.begin(), vec.end()));
            break;
        }
        default:break;                                                  
    }
    if (kIndex==-1) subDim = totDim;
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

void Basis::saveBasis(std::string basisfile) {
    std::ofstream outfile;
    ind_int *d_pt = indexList.data();
    ind_int size = indexList.size();
    save<ind_int>(d_pt, size, &outfile, basisfile);
}

void Basis::saveBasis(std::string basisfile, std::string normfile) {
    std::ofstream outfile;
    save<ind_int>(indexList.data(), (ind_int)indexList.size(), &outfile, basisfile);
    save<double>(normList.data(), (ind_int)normList.size(), &outfile, normfile);
}

ind_int Basis::vecToInd(VecI& v) const {
    ind_int result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < N; i++){if (v.at(i) != 0) result += v[i] * mul[i];}
    return result;
}

pairIndex Basis::vecToInd(VecI& v, VecI& vp) const {
    ind_int r=0, rp=0;
    #pragma omp parallel for reduction(+:r,rp)
    for (int i = 0; i < N; i++){
        if (v.at(i) != 0) r += mul[i];
        if (vp.at(i) != 0) rp += mul[i];
    }
    pairIndex pairInd(r,rp);
    return pairInd;
}

void Basis::indToVec(ind_int index, VecI& v) const {
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

void Basis::indToVec(pairIndex pairInd, VecI& v, VecI& vp) const {
    indToVec(pairInd.first,v);
    indToVec(pairInd.second,vp);
}

void Basis::indToVec(ind_int index, VecI& v, VecI& vp) const {
    indToVec(std::make_pair(fIndexList.at(index/sDim), sIndexList.at(index%sDim)), v, vp);
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
bool Basis::search(ind_int index, ind_int &ind, const std::vector<ind_int> &indList) const {
    // Binary search
    ind_int list_size = indList.size();
    ind_int begin = 0, middle, end = list_size;
    while (end >= begin){
        middle = (end + begin) / 2;
        // in case index>indList[subdim-1], make sure middle<=subdim-1
        if (middle>=list_size) break;
        if (index > indList.at(middle)){
            begin = middle + 1;
            continue;
        }
        else if (index < indList.at(middle)){
            // in case we are using unsigned integer and index<indList[0], make sure middle-1 is not negative
            if (middle>=1){
                end = middle - 1;
                continue;
            }else{
                break;
            }
        }
        else if (index == indList.at(middle)){
            ind = middle;
            return true;
        }
    }
    if (kIndex==-1){
        std::cerr<<"index:"<<index<<" not found in indList (full hilbert space kIndex=-1)!"<<std::endl;
        exit(EXIT_FAILURE);
    }
    return false;
}

bool Basis::search(ind_int index, ind_int &ind) const {
    return search(index, ind, indexList);
}

bool Basis::search(pairIndex pairInd, ind_int &ind) const {
    assert(model != HEISENBERG);
    ind_int fInd, sInd;
    if (model==HUBBARD and kIndex==-1){
        if((search(pairInd.first, fInd, fIndexList)) and (search(pairInd.second, sInd, sIndexList))){
            ind = fInd * sDim + sInd;
            return true;
        }
    }else{
        if((search(pairInd.first, fInd, fIndexList)) and (search(pairInd.second, sInd, sIndexList))){
            ind_int initInd = fInd * sDim + sInd;
            return search(initInd, ind, indexList);
        }
    }
    return false;
}

ind_int Basis::search(ind_int index, const std::vector<ind_int> &indList) const {
    // Binary search
    ind_int list_size = indList.size();
    ind_int begin = 0, middle, end = list_size;
    while (end >= begin){
        middle = (end + begin) / 2;
        if (middle>=list_size) break;
        if (index > indList.at(middle)){
            begin = middle + 1;
            continue;
        }
        else if (index < indList.at(middle)){
            if (middle>=1){
                end = middle -1;
                continue;
            }else{
                break;
            }
        }
        else if (index == indList.at(middle)){
            return middle;
        }
    }
    std::cerr<<"index not found in indList!"<<std::endl;
    exit(EXIT_FAILURE);
}
ind_int Basis::search(ind_int index) const {
    return search(index, indexList);
}
ind_int Basis::search(pairIndex pairInd) const {
    ind_int fInd = search(pairInd.first, fIndexList);
    ind_int sInd = search(pairInd.second, sIndexList);
    return (fInd * sDim + sInd);
}

/*
    *******************************************
    * Implementation for translation symmetry *
    *******************************************
*/
bool Basis::isRep(ind_int initInd, double& norm) const {
    // full hilbert space
    if (kIndex==-1) {
        norm = 1.0;
        return true;
    }
    // smallest index in the cycle?
    switch(model){
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::t_J:{
            pairIndex pairRepI;
            getRepI(initInd,pairRepI);
            // project out double occp
            if (model==LATTICE_MODEL::t_J and (pairRepI.first & pairRepI.second)){
                return false;
            }
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    pairIndex pairRepIp{0,0};
                    for (int i = 0; i < pt_lattice->getOrbNum(); i++){
                        if(bitTest(pairRepI.first,i))bitSet(pairRepIp.first,pt_lattice->getOrbTran(r,i));
                    }
                    if (pairRepI.first > pairRepIp.first) return false;
                    else if (pairRepI.first == pairRepIp.first){
                        for (int j = 0; j < pt_lattice->getOrbNum(); j++){
                            if(bitTest(pairRepI.second,j))bitSet(pairRepIp.second,pt_lattice->getOrbTran(r,j));                            
                        }
                        if (pairRepI.second > pairRepIp.second) return false;
                    }
                }
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                        pairIndex pairRepIp{0,0};
                        for (int i = 0; i < pt_lattice->getOrbNum();i++){
                            if(bitTest(pairRepI.first,i))bitSet(pairRepIp.first,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                        }
                        if (pairRepI.first > pairRepIp.first) return false;
                        else if (pairRepI.first == pairRepIp.first){
                            for (int i = 0; i < pt_lattice->getOrbNum(); i++){
                                if(bitTest(pairRepI.second,i))bitSet(pairRepIp.second,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));                      
                            }
                            if (pairRepI.second > pairRepIp.second) return false;
                        }
                    }    
                }
            }
            break;
        }
        case LATTICE_MODEL::HEISENBERG:{
            if(siteDim>2){
                VecI initVec(N), finalVec(N), finalVec1(N);
                indToVec(initInd, initVec);
                if(PGRepIndex==-1){
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        if (initInd > vecToInd(finalVec)) return false;
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
                            if (initInd > vecToInd(finalVec1)) return false;
                        }        
                    }
                }
            }
            else if(siteDim==2){
                ind_int repIp{0};
                if(PGRepIndex==-1){
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            if(bitTest(initInd,i))bitSet(repIp,pt_lattice->getOrbTran(r,i));
                        }
                        if (initInd > repIp) return false;
                    }
                }else{
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex); p++){
                            for (int i = 0; i < pt_lattice->getOrbNum();i++){
                                if(bitTest(initInd,i)) bitSet(repIp,pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                            }
                            if (initInd > repIp) return false;
                        }        
                    }
                }
            }
            break;
        }
    }
    // norm > 0?
    norm = Norm(initInd);
    // sqrt(infinitismal)>>infinitesimal
    if (std::real(norm)*std::real(norm)>INFINITESIMAL){
        return true;
    }
    return false;
}

double Basis::Norm(ind_int initInd) const {
    if (kIndex==-1) return 1.0;
    cdouble norm = 0.0;
    switch(model){
        case HUBBARD:case t_J:{
            VecI initVec(N), initVecp(N), finalVec(N), finalVecp(N);
            pairIndex pairInd = std::make_pair(fIndexList.at(initInd/sDim), sIndexList.at(initInd%sDim));
            indToVec(pairInd, initVec, initVecp);
            VecI seq, seqp;
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    seq.clear();
                    seqp.clear();
                    for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                        finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        finalVecp.at(pt_lattice->getOrbTran(r,i)) = initVecp.at(i);
                        if (initVec.at(i)==1) seq.push_back(pt_lattice->getOrbTran(r,i));
                        if (initVecp.at(i)==1) seqp.push_back(pt_lattice->getOrbTran(r,i));
                    }
                    if ((vecToInd(finalVec)==pairInd.first) && (vecToInd(finalVecp)==pairInd.second)) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        seq.clear();
                        seqp.clear();
                        for(int i = 0; i < pt_lattice->getOrbNum(); i++){
                            finalVec.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVec.at(i);
                            finalVecp.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVecp.at(i);
                            if (initVec.at(i)==1) seq.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                            if (initVecp.at(i)==1) seqp.push_back(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i)));
                        } 
                        if ((vecToInd(finalVec)==pairInd.first) && (vecToInd(finalVecp)==pairInd.second)) norm += seqSign(seq) * seqSign(seqp) * pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }
                }
                norm /= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
            break;
        }

        case HEISENBERG:{
            VecI initVec(N), finalVec(N);
            indToVec(initInd, initVec);
            if(PGRepIndex==-1){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int i = 0; i < N; i++){
                        finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                    }
                    if (vecToInd(finalVec)==initInd) norm += pt_lattice->expKR(kIndex,r);
                }
                norm /= pt_lattice->getSiteNum();
            }else{
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                        for(int i = 0; i < N; i++){
                            finalVec.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVec.at(i);
                        }
                        if (vecToInd(finalVec)==initInd) norm += pt_lattice->expKR(kIndex,r) * pt_lattice->getChi(PGRepIndex,p);
                    }   
                }
                norm /= pt_lattice->getSiteNum()* pt_lattice->getPGOpNum(PGRepIndex);
            }
            assert(std::abs(std::imag(norm))<INFINITESIMAL);
            return std::sqrt(std::real(norm));
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
            getRepI(ind,repI);
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
            if(siteDim==2){
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    ind_int indp{0};
                    for(int i = 0; i<N; i++){
                        if(bitTest(ind,i))bitSet(indp,pt_lattice->getOrbTran(r,i));
                    }
                    finalInd.push_back(indp);
                }
            }
            else{
                VecI veci(N), vecf(N);
                indToVec(ind, veci);
                for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                    for(int i = 0; i<N; i++){
                        vecf.at(pt_lattice->getOrbTran(r,i)) = veci.at(i);
                    }
                    finalInd.push_back(vecToInd(vecf));
                }
            }
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
        case LATTICE_MODEL::HUBBARD:case LATTICE_MODEL::t_J:{
            pairIndex pairRepI;
            getRepI(repI,pairRepI);
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

        case LATTICE_MODEL::HEISENBERG:{
            if(siteDim==2){
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
            }
            else{
                VecI initVec(N), finalVec(N);
                indToVec(repI, initVec);
                if(PGRepIndex==-1){
                    initNorm *= pt_lattice->getSiteNum();
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for(int i = 0; i < N; i++){
                            finalVec.at(pt_lattice->getOrbTran(r,i)) = initVec.at(i);
                        }
                        finalInd.push_back(vecToInd(finalVec));
                        factorList.push_back(pt_lattice->expKR(kIndex,r)/initNorm);
                    }
                }else{
                    initNorm *= pt_lattice->getSiteNum() * pt_lattice->getPGOpNum(PGRepIndex);
                    for (int r = 0; r < pt_lattice->getSiteNum(); r++){
                        for (int p = 0; p < pt_lattice->getPGOpNum(PGRepIndex);p++){
                            for(int i = 0; i < N; i++){
                                finalVec.at(pt_lattice->getOrbPG(p,pt_lattice->getOrbTran(r,i))) = initVec.at(i);
                            }
                            finalInd.push_back(vecToInd(finalVec));
                            factorList.push_back(pt_lattice->expKR(kIndex,r)* pt_lattice->getChi(PGRepIndex,p)/initNorm);
                        }
                    }
                }
            }

            break;
        }
        
        default:{
            std::cout<<"model not defined! must be: HUBBARD,t_J,HEISENBERG."<<std::endl;
            exit(1);
            break;
        }
    }
}