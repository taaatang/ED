//
//  Operators.cpp
//  ED
//
//  Created by tatang on 10/26/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "Operators.hpp"

/*
    *********************
    * Hamiltonian Class *
    *********************
*/
void Current::row(idx_t rowID, std::vector<MAP<cdouble>>& rowMaps){
    // off diagonal part
    int idx = 0; if(plz=="x") idx = 0; else if(plz=="y") idx = 1; else if(plz=="z") idx=2; else exit(1);
    std::vector<pairIndex> finalIndList;
    std::vector<cdouble> factorList;
    pt_Basis->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < finalIndList.size(); ++i){
        auto pairRepI = finalIndList[i];
        bool isfminRep = pt_Basis->isfMin(pairRepI.first);
        for (auto linkit = Links.begin(); linkit != Links.end(); linkit++){
            int matID = (*linkit).getmatid();
            cdouble factor = CPLX_I * factorList.at(i) * (*linkit).getVal();
            for (auto bondit = (*linkit).begin(); bondit != (*linkit).end(); bondit++){
                auto bondid = bondit - (*linkit).begin();
                auto linkVec = linkit->getvec(linkit->getvecid(bondid));
                double sign = -1.0; if(linkVec.at(idx)<0.) sign = 1.0; else if(linkVec.at(idx)==0.) sign = 0.0;
                cdouble factor_s = factor*sign;
                int siteI = (*bondit).at(0);
                int siteJ = (*bondit).at(1);
                // cp.siteI * cm.siteJ
                cpcm(SPIN::UP, siteI, siteJ, factor_s, finalIndList[i], &rowMaps[matID]);
                cpcm(SPIN::UP, siteJ, siteI, -factor_s, finalIndList[i], &rowMaps[matID]);
                if(isfminRep){
                    cpcm(SPIN::DOWN, siteI, siteJ, factor_s, finalIndList[i], &rowMaps[matID]);
                    cpcm(SPIN::DOWN, siteJ, siteI, -factor_s, finalIndList[i], &rowMaps[matID]);   
                } 
            }
        }
    }
}

Nocc::Nocc(Geometry *pt_lat, Basis *pt_Ba):\
pt_lattice(pt_lat),pt_Basis(pt_Ba),\
SparseMatrix<double>(pt_Ba, pt_Ba, pt_Ba->getSubDim(), 0, pt_lat->getUnitOrbNum()) {
    for(int i = 0; i < pt_lat->getUnitOrbNum(); ++i) diagValList.resize(BaseMatrix<double>::nloc);
    records.resize(pt_lat->getUnitOrbNum());
}

void Nocc::row(idx_t rowID) {
    // diagonal part. occupancy
    VecI occ;
    idx_t repI = pt_Basis->getRepI(rowID);
    pairIndex pairRepI = pt_Basis->getPairRepI(repI);
    pt_lattice->orbOCC(pairRepI, occ);
    idx_t loc_rowID = rowID - BaseMatrix<double>::startRow;
    for (int i = 0; i < pt_lattice->getUnitOrbNum(); ++i) {SparseMatrix<double>::diagValList[i][loc_rowID] = occ[i];}
}

void Nocc::genMat() {
    #pragma omp parallel for
    for(idx_t rowID = BaseMatrix<double>::startRow; rowID < BaseMatrix<double>::endRow; rowID++) row(rowID);
}

double Nocc::count(ORBITAL orbital, dataType* vec) {
    VecI ids = pt_lattice->getOrbID(orbital);
    double sum_final = 0.0;
    for (auto id : ids) {
        double sum = 0.0, part_sum = 0.0;
        #pragma omp parallel for reduction(+:part_sum)
        for(idx_t i = 0; i < BaseMatrix<double>::nloc; ++i){
            part_sum += diagValList[id][i] * (std::real(vec[i]) * std::real(vec[i]) + std::imag(vec[i]) * std::imag(vec[i]));
        }
        MPI_Allreduce(&part_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sum_final += sum;
    }
    return sum_final;
}

double Nocc::count(int id, dataType* vec) {
    double sum = 0.0, part_sum = 0.0;
    #pragma omp parallel for reduction(+:part_sum)
    for(idx_t i = 0; i < BaseMatrix<double>::nloc; ++i){
        part_sum += diagValList[id][i] * (std::real(vec[i]) * std::real(vec[i]) + std::imag(vec[i]) * std::imag(vec[i]));
    }
    MPI_Allreduce(&part_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

void Nocc::count(dataType* vec) {
    for (int id = 0; id < pt_lattice->getUnitOrbNum(); ++id) {
        records.at(id).push_back(count(id, vec));
    }
}

void Nocc::clear( ) {
    records.clear();
}

void Nocc::save(std::string dir) {
    for (int id = 0; id < pt_lattice->getUnitOrbNum(); ++id) {
        ::save(records.at(id).data(), (int)records.at(id).size(), dir + "/orb" + tostr(id));
    }
}
