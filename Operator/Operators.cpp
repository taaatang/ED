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

Current::Current(Geometry *latt, Basis *Bi, Basis *Bf, bool commuteWithSymm):OperatorBase<cdouble>(latt, Bi, Bf, commuteWithSymm) {
    
}

void Current::setDirection(const std::string& plz)
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

void Current::setDirection(Vec3d d) {
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

void Current::print(std::ostream &os) const {
    for (const auto& link : this->myHoppingT) {
        link.print(false, os);
    }
}

void Current::row(idx_t rowID, std::vector<MAP<cdouble>>& rowMaps) {
    std::vector<pairIdx_t> finalIndList;
    std::vector<cdouble> factorList;
    Bf->genSymm(rowID, finalIndList, factorList);
    for (int i = 0; i < (int)finalIndList.size(); ++i) {
        auto pairRepI = finalIndList[i];
        bool isfminRep = Bi->isfMin(pairRepI.first);
        for (const auto& link : this->myHoppingT) {
            int matID = link.getmatid();
            cdouble factor = CPLX_I * factorList.at(i) * link.getVal();
            for (const auto& bond : link.bond()){
                int siteI = bond.at(0);
                int siteJ = bond.at(1);
                // cp.siteI * cm.siteJ
                cpcm(SPIN::UP, siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
                cpcm(SPIN::UP, siteJ, siteI, -factor, pairRepI, &rowMaps[matID]);
                if(isfminRep){
                    cpcm(SPIN::DOWN, siteI, siteJ, factor, pairRepI, &rowMaps[matID]);
                    cpcm(SPIN::DOWN, siteJ, siteI, -factor, pairRepI, &rowMaps[matID]);   
                } 
            }
        }
    }
}

Nocc::Nocc(Geometry *latt, Basis *Bi, Basis *Bf, bool commuteWithSymm):OperatorBase<double>(latt, Bi, Bf, commuteWithSymm, 0, latt->getUnitOrbNum()) {
    for(int i = 0; i < latt->getUnitOrbNum(); ++i) diagValList.resize(BaseMatrix<double>::nloc);
    records.resize(latt->getUnitOrbNum());
}

void Nocc::row(idx_t rowID) {
    // diagonal part. occupancy
    VecI occ;
    idx_t repI = Bf->getRepI(rowID);
    pairIdx_t pairRepI = Bf->getPairRepI(repI);
    latt->orbOCC(pairRepI, occ);
    idx_t loc_rowID = rowID - BaseMatrix<double>::startRow;
    for (int i = 0; i < latt->getUnitOrbNum(); ++i) {
        SparseMatrix<double>::diagValList[i][loc_rowID] = occ[i];
    }
}

void Nocc::construct( ) {
    #pragma omp parallel for
    for(idx_t rowID = BaseMatrix<double>::startRow; rowID < BaseMatrix<double>::endRow; rowID++) {
        row(rowID);
    }
}

double Nocc::count(ORBITAL orbital, dataType* vec) {
    VecI ids = latt->getOrbID(orbital);
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
    double sum = 0.0;
    double part_sum = 0.0;
    #pragma omp parallel for reduction(+:part_sum)
    for(idx_t i = 0; i < BaseMatrix<double>::nloc; ++i){
        part_sum += diagValList[id][i] * (std::real(vec[i]) * std::real(vec[i]) + std::imag(vec[i]) * std::imag(vec[i]));
    }
    MPI_Allreduce(&part_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

void Nocc::count(dataType* vec) {
    for (int id = 0; id < latt->getUnitOrbNum(); ++id) {
        records.at(id).push_back(count(id, vec));
    }
}

void Nocc::clear( ) {
    records.clear();
}

void Nocc::save(const std::string &dir) {
    for (int id = 0; id < latt->getUnitOrbNum(); ++id) {
        ::save(records.at(id).data(), (int)records.at(id).size(), dir + "/orb" + tostr(id));
    }
}
