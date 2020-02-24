//
//  SpectraMain.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//


#include "globalPara.hpp"
#include "Geometry.hpp"
#include "Basis.hpp"
#include "Operators.hpp"
#include "utils.hpp"
#include "PARPACKSolver.hpp"
#include "LANCZOSIterator.hpp"
#include "Spectra.hpp"

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
#include <stdlib.h> // system
#include <chrono>
#include <mpi.h>


int main(int argc, const char * argv[]) {
/*
    ****************************
    * Input And Initialization *
    ****************************
*/
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    OMP_Info(workerID);
    bool COMPUTE_SS = true;
    bool COMPUTE_SQW = true;
    a_int nev = 1;
    // int Nx = 6, Ny = 6;
    // int N = Nx * Ny;
    int N = 36;
    double J1 = 1.0;
    double J2 = 0.0, dJ2 = 0.01;
    int kIndex = 0; // Gamma Point
    std::ifstream infile("spectra_input.txt");
    infile>>N>>kIndex;
    infile.close(); 

    Timer timer;

    bool BASIS_IS_SAVED = false;
    // data directory
    // std::string subDir = std::to_string(Nx) + "by" + std::to_string(Ny);
    std::string subDir = std::to_string(N);
    std::string basisDir = PROJECT_DATA_PATH+"/" + subDir + "/kSpace/Basis/"+std::to_string(kIndex);
    std::string basisfile = basisDir + "/basis";
    std::string normfile = basisDir + "/norm";
    /*
        **************************
        * Lattice Initialization *
        **************************
    */
    TriAngLattice Lattice(N);
    Lattice.addOrb({}); 
    Lattice.construct();
    /*
        **********************
        * Basis Construction *
        **********************
    */
    timer.tik();
    int siteDim = 2;
    VecI occList{N/2, N-N/2};
    Basis B(LATTICE_MODEL::HEISENBERG, &Lattice, occList, kIndex);
    std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED) B.gen(basisfile, normfile);
    else B.gen();
    timer.tok();
    if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kInd ="<<kIndex<<", size="<<B.getSubDim()<<"/"<<B.getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl;
    if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
    /*
        ****************************
        * Hamiltonian Construction *
        ****************************
    */
    timer.tik();
    Link<dataType> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0);
    Link<dataType> J2Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, false);
    J1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
    J2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
    Heisenberg<dataType> H(&Lattice, &B, 2);
    H.pushLinks({J1Link}).pushLinks({J2Link});
    H.genMatPara();
    timer.tok();
    std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.get_nloc()<<"/"<<H.get_dim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
            <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

for(int J2_num = 0; J2_num<1; J2_num++){
    J2 = dJ2 * J2_num;
    std::ofstream outfile;
    // data directory
    // std::string subDir = std::to_string(Nx) + "by" + std::to_string(Ny);
    std::string subDir = std::to_string(N);
    std::string dataDir = PROJECT_DATA_PATH+"/"+ subDir +"/kSpace/Spectra/J2_"+std::to_string(J2_num);
    if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
    if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;
    H.setVal(J2Link.getmatid(),J2);
    /*
        *********************************
        * Diagonalization Using PARPACK *
        *********************************
    */
    PARPACKComplexSolver<double> PDiag(&H, nev);
    PDiag.diag();
    // Ground State
    cdouble* w0 = PDiag.getEigval();
    dataType* gstate = PDiag.getEigvec();
    MPI_Barrier(MPI_COMM_WORLD);
/*
    *************************************
    * Static Spin-Spin Structure Factor *
    *************************************
*/
    if (COMPUTE_SS){
        timer.tik();
        SSOp<dataType> SS(&Lattice,&B);
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin SS ..."<<std::endl<<"********************"<<std::endl;
        std::vector<dataType> vecTmp(SS.get_nlocmax());
        cdouble val;
        std::vector<cdouble> ssvals;
        for (int i = 0; i < Lattice.getSiteNum(); i++){
            val = 0.0;
            SS.setr(i);
            SS.genMatPara();
            SS.MxV(gstate, vecTmp.data());
            vConjDotv<dataType, dataType>(gstate, vecTmp.data(), &val, SS.get_nloc());
            val /= Lattice.getSiteNum();
            ssvals.push_back(val);
            if (workerID==MPI_MASTER) std::cout<<"SS "<<i<<" finished:"<<val<<std::endl;
        }
        // save ss(i)
        if (workerID==MPI_MASTER) save<cdouble>(ssvals.data(), Lattice.getSiteNum(), &outfile, dataDir + "/spinspin");
        timer.tok();
        if (workerID==MPI_MASTER) std::cout<<"SS time:"<<timer.elapse()<<" milliseconds."<<std::endl;
    }
    

/*
    **************************************
    * Dynamic Spin-Spin Structure Factor *
    **************************************
*/
    if (COMPUTE_SQW){
        int krylovDim = 250;
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Sqw ..."<<std::endl<<"********************"<<std::endl;
        for (int kIndexf = 0; kIndexf < Lattice.getSiteNum(); kIndexf++){
            timer.tik();
            if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<kIndexf<<std::endl<<"********************"<<std::endl;
            std::string basisDirp = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/Basis/"+std::to_string(kIndexf);
            std::string basisfilep = basisDirp + "/basis";
            std::string normfilep = basisDirp + "/norm";
            Basis Bp = Basis(LATTICE_MODEL::HEISENBERG, &Lattice, occList, kIndexf);
            if (BASIS_IS_SAVED) {Bp.gen(basisfilep, normfilep);
            }else{Bp.gen();}
            // <Bp|Szq|B>, q = k_B - k_Bp
            SzkOp<dataType> Szq(&Lattice, &B, &Bp);
            Szq.genMatPara();
            Heisenberg<dataType> Hp(&Lattice, &Bp, 2);
            Hp.pushLink(J1Link).pushLink(J2Link);
            Hp.genMatPara();
            Hp.setVal(J2Link, J2);
            timer.tok();
            if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<Hp.get_nloc()<<"/"<<Hp.get_dim()<<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
            SPECTRASolver<dataType> spectra(&Hp, w0[0], &Szq, gstate, H.get_dim(), krylovDim);
            spectra.compute();
            // save alpha, beta
            if (workerID==MPI_MASTER){
                std::string dataPath = dataDir + "/k" + std::to_string(kIndex) + "_kp" + std::to_string(kIndexf);
                system(("mkdir -p " + dataPath).c_str());
                spectra.saveData(dataPath);
            } 
            timer.tok();
            if (workerID==MPI_MASTER) std::cout<<"Sqw time:"<<timer.elapse()<<" milliseconds."<<std::endl<<std::endl;
        }
    }
}   
/*
    *******
    * End *
    *******
*/
    MPI_Finalize();
    return 0;
}

