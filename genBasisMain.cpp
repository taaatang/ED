//
//  genBasisMain.cpp
//  ED
//
//  Created by tatang on 12/17/19.
//  Copyright © 2019 tatang. All rights reserved.
//
#include "globalPara.hpp"
#include "Geometry.hpp"
#include "Basis.hpp"
#include "utils.hpp"
#include "HelperClass.hpp"

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
#include <stdlib.h> // system
#include <chrono>
#include <mpi.h>
#ifdef OMP_
    #include <omp.h>
#endif

int main(int argc, const char * argv[]) {
/*
    ********************
    * MPI and OMP Info *
    ********************
*/
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    
    OMP_Info(workerID);
    
/*
    ****************************
    * Input And Initialization *
    ****************************
*/
    int Nx = 4, Ny = 4;
    int N = Nx * Ny;
    int kIndex = 0;
    Timer timer;
    
    if (workerNum > 1){
        kIndex = workerID;
    }

    if (kIndex < N){
        // create data path
        std::string dataDir = PROJECT_DATA_PATH+"/"+std::to_string(Nx)+"by"+std::to_string(Ny)+"/kSpace/Basis/"+std::to_string(kIndex);
        system(("mkdir -p " + dataDir).c_str());
        std::string basisfile = dataDir + "/basis";
        std::string normfile = dataDir + "/norm";
        
    /*
        ************************************
        * Lattice and Basis Initialization *
        ************************************
    */
        // geometry class
        TriAngLattice Lattice(Nx,Ny);
        Lattice.addOrb({});
        Lattice.construct();
        int siteDim = 2;
        VecI occList{N/2, N/2};
    
        timer.tik();
        Basis B(LATTICE_MODEL::HEISENBERG, &Lattice, occList, kIndex);
        B.gen();
        timer.tok();
        std::cout<<"WorkerID:"<<workerID<<", kInd ="<<B.getkIndex()<<", size="<<B.getSubDim()<<"/"<<B.getTotDim()<<". Basis construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
        B.saveBasis(basisfile, normfile);
        // check
        // std::cout<<"WorkerID"<<workerID<<", begin reading basis..."<<std::endl;
        // tic = std::chrono::system_clock::now();
        // Basis Bp(MODEL,siteDim, dimList, &Lattice);
        // Bp.gen(kIndex, basisfile, normfile);
        // toc = std::chrono::system_clock::now();
        // elapsed_seconds = toc - tic;

        // std::cout<<"WorkerID:"<<workerID<<", kInd ="<<Bp.kIndex<<", size="<<Bp.subDim<<"/"<<Bp.totDim<<". Basis construction time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl<<"Bp.kRepList size:"<<Bp.indexList.size()<<", Bp.norm size:"<<Bp.normList.size()<<std::endl;
        // for (ind_int i = 0; i < B.subDim; i++){
        //     if ((Bp.indexList.at(i)!=B.indexList.at(i)) or (Bp.normList.at(i)!=B.normList.at(i))){
        //         std::cout<<"workerID:"<<workerID<<", inconsistent data at i="<<i<<std::endl;
        //         exit(1);
        //     }
        // }
        //  std::cout<<"workerID:"<<workerID<<", B and Bp maches!"<<std::endl;
    }

/*
    *******
    * End *
    *******
*/
    MPI_Finalize();
    return 0;
}



