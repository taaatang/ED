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


// #ifdef SAXPY
//     // generate Hamiltonian in the subspacd labeled by kIndex
//     void Heisenberg::genSubMatMap(int kIndex, int couplingNum, int polarNum){
//         assert(couplingNum==matrix_num);
//         if (kIndex==-1){
//             genMat(pt_lattice, pt_Basis, couplingNum, polarNum);
//             return;
//         }
//         clear();
//         for (int matID = 0; matID < matrix_num; matID++){
//             pt_rowInitList[matID]->reserve(nloc+1);
//             pt_valList[matID]->reserve(nloc*polarNum*pt_lattice->getOrbNum()/2+nloc);
//             pt_colList[matID]->reserve(nloc*polarNum*pt_lattice->getOrbNum()/2+nloc);
//         }
//         MAP rowMap;
//         MAPIT it;

//         std::vector<double> k(2,0.0);
//         pt_lattice->getK(kIndex, k.data());
//         ind_int *finalIndList_ = new(std::nothrow) ind_int[pt_lattice->getOrbNum()]; assert(finalIndList_!=NULL);
//         int *initVec_ = new(std::nothrow) int[pt_lattice->getOrbNum()]; assert(initVec_!=NULL);
//         ind_int counter;
//         double initNorm, finalNorm;
        
//         // calculate <R1k|H*Pk|R2k>/norm1/norm2
//         for (int couplingID = 0; couplingID < couplingNum; couplingID++){
//             counter = 0;
//             pt_rowInitList[couplingID]->push_back(counter);
//             for (ind_int rowID = startRow; rowID < endRow; rowID++){
//                 // clean map
//                 rowMap.clear();
//                 pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList_);
//                 initNorm = pt_Basis->getNorm(rowID);
//                 for (int i = 0; i < pt_lattice->getOrbNum(); i++){
//                     if (DEBUG) std::cout<<"trans vec:"<<i<<" begins..."<<std::endl;
//                     pt_Basis->indToVec(finalIndList_[i], initVec_);
//                     std::vector<double> r(2,0.0);
//                     pt_lattice->getOrbR(i, r.data());
//                     cdouble factor = std::exp(2*PI*CPLX_I*(k[0]*r[0]+k[1]*r[1]))/pt_lattice->getSiteNum()/initNorm;
//                     for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                         for (int polarID = 0; polarID < polarNum; polarID++){
//                             int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
//                             // sz.siteID * sz.siteIDp
//                             szsz(siteID, siteIDp, parameters.at(couplingID)*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
//                             // 1/2 * sm.siteID * sp.siteIDp
//                             spsm(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
//                             // 1/2 * sp.siteID * sm.siteIDp
//                             smsp(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList_[i], initVec_, pt_Basis, &rowMap);
//                         }
//                     }
//                 }
//                 newRow(couplingID, &rowMap, counter);
//             }
//         }
//         delete [] initVec_;
//         delete [] finalIndList_;
//     }

//     void Heisenberg::genMat(int couplingNum, int polarNum){
//         clear();
//         for (int matID = 0; matID < matrix_num; matID++){
//             pt_rowInitList[matID]->reserve(nloc+1);
//             pt_valList[matID]->reserve(nloc*polarNum*pt_lattice->getOrbNum()/2+nloc);
//             pt_colList[matID]->reserve(nloc*polarNum*pt_lattice->getOrbNum()/2+nloc);
//         }
//         MAP rowMap;
//         MAPIT it;

//         ind_int counter;
//         int *initVec_ = new int[pt_lattice->getOrbNum()];
//         for (int couplingID = 0; couplingID < couplingNum; couplingID++){
//             counter = 0;
//             pt_rowInitList[couplingID]->push_back(counter);
//             for (ind_int rowID = startRow; rowID < endRow; rowID++){
//                 rowMap.clear();
//                 pt_Basis->indToVec(pt_Basis->indexList.at(rowID), initVec_);
//                 for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                     for (int polarID = 0; polarID < polarNum; polarID++){
//                         int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
//                         // sz.siteID * sz.siteIDp
//                         szsz(siteID, siteIDp, parameters.at(couplingID), rowID, initVec_, pt_Basis, &rowMap);
//                         // 1/2 * sm.siteID * sp.siteIDp
//                         spsm(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
//                         // 1/2 * sp.siteID * sm.siteIDp
//                         smsp(siteID, siteIDp, parameters.at(couplingID)/2.0, rowID, initVec_, pt_Basis, &rowMap);
//                     }
//                 }
//                 newRow(couplingID, &rowMap, counter);
//             }
//         }
//         delete [] initVec_;
//     }
// #else
//     void Heisenberg::row(int kIndex, int couplingNum, int polarNum, ind_int rowID, MAP *rowMap){
//         rowMap->clear();
//         std::vector<ind_int> finalIndList(pt_lattice->getOrbNum());
//         std::vector<int> initVec(pt_lattice->getOrbNum());
//         double initNorm, finalNorm;
//         // Bug Here! different threads may change the common initVec in the Basis object.
//         pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList.data());
//         double kx = pt_lattice->KLattice_[kIndex].coord[0];
//         double ky = pt_lattice->KLattice_[kIndex].coord[1];
//         initNorm = pt_Basis->getNorm(rowID);
//         for (int i = 0; i < pt_lattice->getOrbNum(); i++){
//             pt_Basis->indToVec(finalIndList[i], initVec.data());
//             int rx = pt_lattice->lattice_[i].coord[0];
//             int ry = pt_lattice->lattice_[i].coord[1];
//             cdouble factor = std::exp(2*PI*CPLX_I*(kx*rx+ky*ry))/pt_lattice->getOrbNum()/initNorm;
//             for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                 for (int couplingID = 0; couplingID < couplingNum; couplingID++){
//                     for (int polarID = 0; polarID < polarNum; polarID++){
//                         int siteIDp = pt_lattice->bondMaps_[couplingID]->at(siteID).at(polarID);
//                         // sz.siteID * sz.siteIDp
//                         szsz(siteID, siteIDp, parameters.at(couplingID)*factor,finalIndList[i], initVec.data(), pt_Basis, rowMap);
//                         // 1/2 * sm.siteID * sp.siteIDp
//                         spsm(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList[i], initVec.data(), pt_Basis, rowMap);
//                         // 1/2 * sp.siteID * sm.siteIDp
//                         smsp(siteID, siteIDp, parameters.at(couplingID)/2.0*factor,finalIndList[i], initVec.data(), pt_Basis, rowMap);
//                     }
//                 }
//             }
//         }
//     }
//     void Heisenberg::pushRow(ind_int start, MAP *rowMap){
//         int count = 0;
//         for (auto it = rowMap->begin(); it != rowMap->end(); it++){
//             colList.at(start+count) = it->first;
//             switch(PARTITION){
//                 case ROW_PARTITION:{
//                     valList.at(start+count) = std::conj(it->second);
//                     break;
//                 }
//                 case COL_PARTITION:{
//                     valList.at(start+count) = it->second;
//                     break;
//                 }
//             }
//             count++;
//         }
//         assert(count==rowMap->size());
//     }
// // generate Hamiltonian in the subspacd labeled by kIndex in parallel
//     void Heisenberg::genSubMatMapPara(int kIndex, int couplingNum, int polarNum){
//         if (kIndex==-1){
//             genMat(couplingNum, polarNum);
//             return;
//         }
//         clear();
//         rowInitList.reserve(nloc+1);
//         colList.reserve(nloc*couplingNum*polarNum*pt_lattice->getOrbNum()/2+nloc);
//         valList.reserve(nloc*couplingNum*polarNum*pt_lattice->getOrbNum()/2+nloc);

//         int threadNum;
//         #pragma omp parallel
//         {
//             #pragma omp master
//             threadNum = omp_get_num_threads();
//         }
        
//         std::vector<MAP*> rowMapList(threadNum);
//         for (int i = 0; i < threadNum; i++) {
//             rowMapList[i] = new MAP;
//             rowMapList[i] ->reserve(couplingNum*polarNum*pt_lattice->getOrbNum()+nloc);
//         }
//         std::vector<ind_int> startList(threadNum);
        
//         ind_int counter = 0;
//         rowInitList.push_back(counter);
        
//         // calculate <R1k|H*Pk|R2k>/norm1/norm2
//         for (ind_int rowID = startRow; rowID < endRow; rowID+=threadNum){
//             #pragma omp parallel shared(rowMapList, startList)
//             {   
//                 assert(omp_get_num_threads()==threadNum);
//                 int threadID = omp_get_thread_num();
//                 ind_int thRowID = rowID + threadID;
//                 bool is_row = thRowID < endRow;
//                 if (is_row) row(kIndex, couplingNum, polarNum, thRowID, rowMapList.at(threadID));
//                 #pragma omp barrier

//                 #pragma omp master
//                 {
//                     for (int i = 0; (i < threadNum)&&(rowID+i<endRow); i++){
//                         startList[i] = counter;
//                         counter += rowMapList.at(i)->size();
//                         rowInitList.push_back(counter);
//                     }
//                     colList.resize(counter);
//                     valList.resize(counter);
//                 }
//                 #pragma omp barrier

//                 if (is_row) pushRow(startList[threadID], rowMapList.at(threadID));
//             }
//         }
//         for (int i = 0; i < threadNum; i++) delete rowMapList[i];
//     }

// generate Hamiltonian in the subspacd labeled by kIndex
void Heisenberg::genMat(){
    int kIndex = pt_Basis->getkIndex();
    if (kIndex==-1){
        genMatFull();
        return;
    }
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
        for (int i = 0; i < pt_lattice->getSiteNum(); i++){
            pt_Basis->indToVec(finalIndList[i], initVec);
            cdouble factor = pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
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
            for (int i = 0; i < pt_lattice->getSiteNum(); i++){
                pt_Basis->indToVec(finalIndList[i], initVec);
                cdouble factor = pt_lattice->expKR(kIndex,i)/pt_lattice->getSiteNum()/initNorm;
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

void Heisenberg::genMatFull(){
    clear();
    // rowInitList.reserve(nloc+1);
    // colList.reserve(nloc*couplingNum*polarNum*pt_lattice->getOrbNum()/2+nloc);
    // valList.reserve(nloc*couplingNum*polarNum*pt_lattice->getOrbNum()/2+nloc);
    MAP rowMap;
    // initialize rowInitList
    for (int i = 0; i < spmNum; i++) pushRow(&rowMap,i);
    dataType dval;
    VecI initVec(pt_lattice->getOrbNum());
    for (ind_int rowID = startRow; rowID < endRow; rowID++){
        rowMap.clear();
        ind_int initInd = pt_Basis->getRepI(rowID);
        pt_Basis->indToVec(initInd, initVec);

        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            double val = (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, val, initInd, initVec, &rowMap);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, val/2.0, initInd, initVec, &rowMap);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, val/2.0, initInd, initVec, &rowMap);
            }
        }
        pushRow(&rowMap);

        for (auto linkit = NCLinks.begin(); linkit != NCLinks.end(); linkit++){
            rowMap.clear();
            double val = (*linkit)->getVal();
            for (auto bondit = (*linkit)->begin(); bondit != (*linkit)->end(); bondit++){
                int siteID = (*bondit).at(0);
                int siteIDp = (*bondit).at(1);
                // sz.siteID * sz.siteIDp
                szsz(siteID, siteIDp, val, initInd, initVec, &rowMap);
                // 1/2 * sm.siteID * sp.siteIDp
                spsm(siteID, siteIDp, val/2.0, initInd, initVec, &rowMap);
                // 1/2 * sp.siteID * sm.siteIDp
                smsp(siteID, siteIDp, val/2.0, initInd, initVec, &rowMap);
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

// SzqOneHalf::SzqOneHalf(ind_int totDim){
//     MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
//     MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
//     dim = totDim;
//     nlocmax = (dim + workerNum - 1)/workerNum;
//     ntot = nlocmax * workerNum;
//     startRow = workerID * nlocmax;
//     endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
//     nloc = endRow - startRow;
// }

// SzqOneHalf::~SzqOneHalf(){}

// void SzqOneHalf::genMat(Geometry* pt_lattice, Basis* pt_Basis, BasisXY q){
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


// // generate
// void SzqOneHalf::genMat(Geometry* pt_lattice, Basis* pt_B1, Basis* pt_B2, BasisXY q){
//     clear();
//     rowInitList.reserve(nloc+1);
//     colList.reserve(nloc);
//     valList.reserve(nloc);
//     cdouble dval;
//     int *initVec_ = new(std::nothrow) int[pt_lattice->getOrbNum()]; assert(initVec_!=NULL);
//     switch(PARTITION){
//         case ROW_PARTITION:{
//             ind_int counter = 0;
//             rowInitList.push_back(counter);
//             ind_int colID;
//             for (ind_int rowID = startRow; rowID < endRow; rowID++){
//                 if (pt_B1->search(pt_B2->indexList.at(rowID),colID)){
//                     dval = 0.0;
//                     pt_B1->indToVec(pt_B1->indexList.at(colID), initVec_);
//                     for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                         dval += szMat.at(initVec_[siteID]) * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
//                     }
//                     dval *= pt_B2->getNorm(rowID)/pt_B1->getNorm(colID);
//                     colList.push_back(colID);
//                     valList.push_back(dval);
//                     counter++;
//                 }
//                 rowInitList.push_back(counter);
//             }
//             break;
//         }
//         case COL_PARTITION:{
//             ind_int counter = 0;
//             rowInitList.push_back(counter);
//             ind_int colID;
//             for (ind_int rowID = startRow; rowID < endRow; rowID++){
//                 if (pt_B1->search(pt_B2->indexList.at(rowID),colID)){
//                     dval = 0.0;
//                     pt_B2->indToVec(pt_B2->indexList.at(rowID), initVec_);
//                     for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                         dval += szMat[initVec_[siteID]] * std::exp(-CPLX_I * (q[0]*pt_lattice->lattice_[siteID].coordxy[0] + q[1]*pt_lattice->lattice_[siteID].coordxy[1]));
//                     }
//                     dval *= pt_B1->getNorm(rowID)/pt_B2->getNorm(colID);
//                     colList.push_back(colID);
//                     valList.push_back(dval);
//                     counter++;
//                 }
//                 rowInitList.push_back(counter);
//             }
//             break;
//         }
//     }
//     delete [] initVec_;
// }

// /*
//     *****************
//     * SS Correlator *
//     *****************
// */
// SSOneHalf::SSOneHalf(ind_int totDim){
//     MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
//     MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
//     dim = totDim;
//     nlocmax = (dim + workerNum - 1)/workerNum;
//     ntot = nlocmax * workerNum;
//     startRow = workerID * nlocmax;
//     endRow = (startRow + nlocmax)<dim?(startRow + nlocmax):dim;
//     nloc = endRow - startRow;
// }

// SSOneHalf::~SSOneHalf(){}

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

// // generate matrix in subsapce labeled by kIndex for sum.r:Sr*Sr+dr, dr is labeled by rIndex
// void SSOneHalf::genPairMat(int kIndex, Geometry* pt_lattice, Basis* pt_Basis, int rIndex){
//     clear();
//     rowInitList.reserve(nloc+1);
//     colList.reserve(nloc*pt_lattice->getOrbNum()/2+nloc);
//     valList.reserve(nloc*pt_lattice->getOrbNum()/2+nloc);
//     MAP rowMap;
//     MAPIT it;

//     double kx = pt_lattice->KLattice_[kIndex].coord[0];
//     double ky = pt_lattice->KLattice_[kIndex].coord[1];
//     int drx = pt_lattice->lattice_[rIndex].coord[0];
//     int dry = pt_lattice->lattice_[rIndex].coord[1];
//     ind_int *finalIndList_ = new(std::nothrow) ind_int[pt_lattice->getOrbNum()]; assert(finalIndList_!=NULL);
//     int *initVec_ = new(std::nothrow) int[pt_lattice->getOrbNum()]; assert(initVec_!=NULL);

//     ind_int counter = 0;
//     cdouble dval;
//     double initNorm, finalNorm;

//     rowInitList.push_back(counter);
//     for (ind_int rowID = startRow; rowID < endRow; rowID++){
//         rowMap.clear();
//         initNorm = pt_Basis->getNorm(rowID);
//         pt_Basis->genTranslation(pt_lattice, pt_Basis->indexList.at(rowID), finalIndList_);
//         for (int i = 0; i < pt_lattice->getOrbNum(); i++){
//             pt_Basis->indToVec(finalIndList_[i], initVec_);
//             dval = 0.0;
//             int rx = pt_lattice->lattice_[i].coord[0];
//             int ry = pt_lattice->lattice_[i].coord[1];
//             cdouble factor = std::exp(2*PI*CPLX_I*(kx*rx+ky*ry))/pt_lattice->getOrbNum()/initNorm;
//             for (int siteID = 0; siteID < pt_lattice->getOrbNum(); siteID++){
//                 int siteIDp;
//                 BasisLat coord;
//                 coord[0] = pt_lattice->lattice_[siteID].coord[0] + drx;
//                 coord[1] = pt_lattice->lattice_[siteID].coord[1] + dry;
//                 assert(pt_lattice->coordToInd(coord, siteIDp));
//                 // sz.siteID * sz.siteIDp
//                 szsz(siteID, siteIDp, factor, finalIndList_[i], initVec_, pt_Basis, &rowMap);
//                 // 1/2 * sm.siteID * sp.siteIDp
//                 spsm(siteID, siteIDp, factor/2.0, finalIndList_[i], initVec_, pt_Basis, &rowMap);
//                 // 1/2 * sp.siteID * sm.siteIDp
//                 smsp(siteID, siteIDp, factor/2.0, finalIndList_[i], initVec_, pt_Basis, &rowMap);
//             }
//         }
//         newRow(&rowMap, counter);
//     }
//     delete [] initVec_;
//     delete [] finalIndList_;
// }

