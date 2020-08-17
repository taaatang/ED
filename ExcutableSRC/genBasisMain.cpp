//
//  genBasisMain.cpp
//  ED
//
//  Created by tatang on 12/17/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "../Basis/Basis.hpp"
#include "../Utils/timer.hpp"

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
#include <stdlib.h> // system
#include <mpi.h>
#include <omp.h>

int main(int argc, const char * argv[]) {
/*
    ********************
    * MPI and OMP Info *
    ********************
*/
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
    
/*
    ****************************
    * Input And Initialization *
    ****************************
*/
    int N, Nx, Ny, Nu, Nd;
    int kIndex;
    std::ifstream infile("../Input/symm_input.txt");
    // infile>>Nx>>Ny>>Nu>>Nd;
    // N = Nx * Ny;
    infile>>kIndex;
    infile.close(); 
    N=36, Nu=18, Nd=18;
    
    std::string subDir = "/"+std::to_string(N);
    // std::string subDir = "/"+std::to_string(Nx)+"by"+std::to_string(Ny)+"_"+std::to_string(Nu)+"u"+std::to_string(Nd)+"d";
    // std::string subDir = "/sq"+std::to_string(N)+"_"+std::to_string(Nu)+"u"+std::to_string(Nd)+"d";
    Timer timer;
    
    if (workerNum > 1){
        kIndex = workerID;
    }

    // for(kIndex=0; kIndex<N; kIndex++){
    if(kIndex < N){
        // create data path
        std::string dataDir = PROJECT_DATA_PATH+ subDir +"/kSpace/Basis/"+std::to_string(kIndex);
        system(("mkdir -p " + dataDir).c_str());
        std::string basisfile = dataDir + "/basis";
        std::string normfile = dataDir + "/norm";
        
    /*
        ************************************
        * Lattice and Basis Initialization *
        ************************************
    */
        // SquareLattice Lattice(N);
        // Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
        // Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
        TriAngLattice Lattice(N);
        Lattice.addOrb({});
        Lattice.construct();
        int siteDim = 2;
        VecI occList{Nu, Nd};
    
        timer.tik();
        Basis B(LATTICE_MODEL::HUBBARD, &Lattice, occList, kIndex);
        B.gen();
        timer.tok();
        std::cout<<"WorkerID:"<<workerID<<", kInd ="<<B.getkIndex()<<", size="<<B.getSubDim()<<"/"<<B.getTotDim()<<". Basis construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
        B.saveBasis(basisfile, normfile);

        //Paralel Construction
        // timer.tik();
        // Basis B(LATTICE_MODEL::HUBBARD, &Lattice, occList, kIndex);
        // B.gen(workerID,workerNum);
        // timer.tok();
        // std::cout<<"WorkerID:"<<workerID<<", kInd ="<<B.getkIndex()<<", local size="<<B.getLocDim()<<"/"<<B.getTotDim()<<". Basis construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
        // int msg = 1;
        // bool is_app = true;
        // if(workerID==0){
        //     B.saveBasis(basisfile, normfile);
        //     std::cout<<"workerID:"<<workerID<<", data saved!"<<std::endl;
        //     if(workerNum>1){
        //         MPI_Send(&msg,1,workerID+1);
        //     }
        // }else if((workerID>0) and (workerID<(workerNum-1))){
        //     MPI_Recv(&msg,1,workerID-1);
        //     B.saveBasis(basisfile, normfile, is_app);
        //     std::cout<<"workerID:"<<workerID<<", data saved!"<<std::endl;
        //     MPI_Send(&msg,1,workerID+1);
        // }else if(workerID==(workerNum-1)){
        //     MPI_Recv(&msg,1,workerID-1);
        //     B.saveBasis(basisfile, normfile, is_app);
        //     std::cout<<"workerID:"<<workerID<<", data saved!"<<std::endl;
        // }

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



