//
//  Operators.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Operators.hpp"

/*
    ************************
    * Spin Operators Class *
    ************************
*/
SpinOperator::SpinOperator(Basis* pt_Ba, LATTICE_MODEL mod, int dim):pt_Basis(pt_Ba),model(mod),spinDim(dim){
    double s = (double)(spinDim - 1)/2.0;
    double m = s;
    for (int i = 0; i < spinDim; i++){
        szMat.push_back(m);
        spMat.push_back(std::sqrt(s*(s+1.0)-m*(m+1.0)));
        smMat.push_back(std::sqrt(s*(s+1.0)-m*(m-1.0)));
        m -= 1.0;
    }
}

/*
    *********************
    * Hamiltonian Class *
    *********************
*/
template <class T>
void Heisenberg::row(ind_int rowID, std::vector<MAP>& rowMaps){
    int kIndex = pt_Basis->getkIndex();
    VecI initVec(pt_lattice->getOrbNum());
    double initNorm, finalNorm;
    std::vector<ind_int> finalIndList;
    pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
    initNorm = pt_Basis->getNorm(rowID);
    for (int i = 0; i < finalIndList.size(); i++){
        pt_Basis->indToVec(finalIndList[i], initVec);
        cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            cdouble factor1 = factor * (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, factor1, finalIndList[i], initVec, &rowMaps[0]);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMaps[0]);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMaps[0]);
            }
        }
    }

    for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
            factor *= (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                int matID = (*linkit)->getmatid();
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMaps[matID]);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMaps[matID]);
            }
        }
    }
}
// generate Hamiltonian in the subspacd labeled by kIndex
template <class T>
void Heisenberg::genMat(){
    int kIndex = pt_Basis->getkIndex();
    clear();
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < spmNum; i++) pushRow(&rowMap,i);
    VecI initVec(pt_lattice->getOrbNum());
    double initNorm, finalNorm;
    // calculate <R1k|H*Pk|R2k>/norm1/norm2
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        std::vector<ind_int> finalIndList;
        pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
        initNorm = pt_Basis->getNorm(rowID);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor1 = factor * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor1, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor1/2.0, finalIndList[i], initVec, &rowMap);
                }
            }
        }
        pushRow(&rowMap);

        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->indToVec(finalIndList[i], initVec);
                cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
                factor *= (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // sz.siteID * sz.siteIDp
                    szsz(siteID, siteIDp, factor, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sm.siteID * sp.siteIDp
                    spsm(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                    // 1/2 * sp.siteID * sm.siteIDp
                    smsp(siteID, siteIDp, factor/2.0, finalIndList[i], initVec, &rowMap);
                }
            }
            pushRow(&rowMap, (*linkit)->getmatid());
        }
    }
}


void Hubbard::row(ind_int rowID, MAP* rowMap, int matID){

}

void Hubbard::genMat(){
    int kIndex = pt_Basis->getkIndex();
    clear();
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < spmNum; i++) pushRow(&rowMap,i);
    VecI initVec(pt_lattice->getOrbNum()), initVecp(pt_lattice->getOrbNum());
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        /*
            *****************
            * Constant Part *
            * ***************
        */
        rowMap.clear();
        // diagonal part. occupancy and double-occ
        VecI occ, docc;
        ind_int initInd = pt_Basis->getRepI(rowID);
        pt_Basis->indToVec(initInd, initVec, initVecp);
        pt_lattice->orbOCC(initVec, initVecp, occ, docc);
        double val = diagVal(occ,docc);
        cdouble diag_val = 0.0;
        // off diagonal part
        std::vector<ind_int> finalIndList;
        std::vector<cdouble> factorList;
        pt_Basis->genTranslation(rowID, finalIndList, factorList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
            if(finalIndList[i]==initInd)diag_val += val*factorList[i];
            for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
                cdouble factor = factorList.at(i) * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteI = (*bondit).at(0);
                    int siteJ = (*bondit).at(1);
                    // cp.siteI * cm.siteJ
                    cpcm(SPIN::SPIN_UP, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_UP, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteI, siteJ, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteJ, siteI, factor, finalIndList[i], initVec, initVecp, &rowMap);   
                }
            }
        }
        diag(rowID,diag_val,&rowMap);
        pushRow(&rowMap);

        /*
            ***********************
            * Time Dependent Part *
            * *********************
        */
        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            for (int i = 0; i < finalIndList.size(); i++){
                pt_Basis->indToVec(finalIndList[i], initVec, initVecp);
                cdouble factor = factorList.at(i) * (*linkit)->getVal();
                for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                    int siteID = (*bondit).at(0);
                    int siteIDp = (*bondit).at(1);
                    // cm.siteID * cp.siteIDp
                    cpcm(SPIN::SPIN_UP, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_UP, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteID, siteIDp, factor, finalIndList[i], initVec, initVecp, &rowMap);
                    cpcm(SPIN::SPIN_DOWN, siteIDp, siteID, factor, finalIndList[i], initVec, initVecp, &rowMap); 
                }
            }
            pushRow(&rowMap, (*linkit)->getmatid());
        }
    }
}
/*
    *******************
    * SzSz Correlator *
    *******************
*/


// void SzkOp::genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q){
//     clear();
//     diagValList.reserve(nloc);
//     cdouble dval;
//     int *initVec_ = new(std::nothrow) int[pt_lattice->getOrbNum()]; assert(initVec_!=NULL);
//     for (ind_int rowID = startRow; rowID < endRow; rowID++){
//         dval = 0.0;
//         pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
//         for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//             dval += szMat.at(initVec_[siteID]) * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
//         }
//         diagValList.push_back(dval);
//     }
//     delete [] initVec_;
// }

/*
    *********
    * Sz(k) *
    *********
*/
// generate
void SzkOp::genMat(){
    clear();
    reserve(1);
    MAP rowMap;
    pushRow(&rowMap);
    // factor[n] =  exp(-i*q*Rn) = exp(i*(Kf-Ki)*Rn)
    std::vector<cdouble> factor(pt_lattice->getOrbNum());
    for (int i = 0; i < pt_lattice->getOrbNum(); i++) factor[i] = pt_lattice->expKR(Kf,i)/pt_lattice->expKR(Ki,i);
    cdouble dval;
    VecI initVec(pt_lattice->getOrbNum());
    switch(PARTITION){
        case ROW_PARTITION:{
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                rowMap.clear();
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bi->indToVec(pt_Bi->getRepI(colID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += getSz(siteID,initVec) * factor[siteID];
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
            ind_int colID;
            for (ind_int rowID = startRow; rowID < endRow; rowID++){
                if (pt_Bi->search(pt_Bf->getRepI(rowID),colID)){
                    dval = 0.0;
                    pt_Bf->indToVec(pt_Bf->getRepI(rowID), initVec);
                    for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
                        dval += getSz(siteID,initVec) * factor[siteID];
                    }
                    dval *= pt_Bi->getNorm(rowID)/pt_Bf->getNorm(colID);
                    rowMap[colID] = dval;
                }
                pushRow(&rowMap);
            }
            break;
        }
    }
}

// /*
//     *****************
//     * SS Correlator *
//     *****************
// */
// generate matrix in subsapce labeled by kIndex for sum.r:Sr*Sr+dr, dr is labeled by rIndex
void SSOp::genPairMat(int rIndex){
    int kIndex = pt_Basis->getkIndex();
    VecI siteJList(pt_lattice->getOrbNum());
    VecD coordi(3), coordr(3), coordf(3);
    pt_lattice->getSiteR(rIndex, coordr.data());
    for (int i = 0; i < pt_lattice->getOrbNum(); i++){
        pt_lattice->getOrbR(i,coordi.data());
        vecXAdd(1.0, coordi.data(), 1.0, coordr.data(), coordf.data(), 3);
        int siteJ;
        if (pt_lattice->coordToOrbid(coordf.data(), siteJ)){
            siteJList[i] = siteJ;
        }else{
            std::cout<<"translation position not found for orbid = "<<i<<", transVecid = "<<rIndex<<std::endl;
            exit(1);
        }
    }
    clear();
    reserve(pt_lattice->getOrbNum()/2+1);
    MAP rowMap;
    pushRow(&rowMap);
    VecI initVec(pt_lattice->getOrbNum());
    double initNorm, finalNorm;
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        initNorm = pt_Basis->getNorm(rowID);
        std::vector<ind_int> finalIndList;
        pt_Basis->genTranslation(pt_Basis->getRepI(rowID), finalIndList);
        for (int i = 0; i < finalIndList.size(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = (kIndex==-1)?1.0:pt_lattice->expKR(pt_Basis->getkIndex(),i)/pt_lattice->getSiteNum()/initNorm;
            for (int siteI = 0; siteI < pt_lattice->getOrbNum(); siteI++){
                int siteJ = siteJList[siteI];
                // sz.siteI * sz.siteJ
                szsz(siteI, siteJ, factor, finalIndList[i], initVec, &rowMap);
                // 1/2 * sm.siteI * sp.siteJ
                spsm(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
                // 1/2 * sp.siteI * sm.siteJ
                smsp(siteI, siteJ, factor/2.0, finalIndList[i], initVec, &rowMap);
            }
        }
        pushRow(&rowMap);
    }
}

// // generate Si*Sj in full hilbert space
// void SSOneHalf::genPairMat(Geometry* pt_lattice, Basis* pt_Basis, int siteID, int initSiteID){
//     clear();
//     rowInitList.reserve(nloc+1);
//     colList.reserve(nloc*3);
//     valList.reserve(nloc*3);
//     MAP rowMap;
//     MAPIT it;

//     ind_int counter = 0;
//     rowInitList.push_back(counter);
//     int *initVec_ = new int[pt_lattice->getOrbNum()];
//     for (ind_int rowID = startRow; rowID < endRow; rowID++){
//         rowMap.clear();
//         pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
//         // sz.siteID * sz.siteIDp
//         szsz(siteID, initSiteID, 1.0, rowID, initVec_, pt_Basis, &rowMap);
//         // 1/2 * sm.siteID * sp.siteIDp
//         spsm(siteID, initSiteID, 0.5, rowID, initVec_, pt_Basis, &rowMap);
//         // 1/2 * sp.siteID * sm.siteIDp
//         smsp(siteID, initSiteID, 0.5, rowID, initVec_, pt_Basis, &rowMap);

//         newRow(&rowMap, counter);
//     }
//     delete [] initVec_;
// }
