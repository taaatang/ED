//
//  main.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//
#define TEST
// #undef TEST
#define OMP_
// #undef OMP_

#include "constant.hpp"
#include "parpackSolver.hpp"
#include "lanczos.hpp"

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
    ****************************
    * Input And Initialization *
    ****************************
*/
    MPI_Init(NULL, NULL);
    int ompThreadsNum;
    a_int nev = 1;
    int Nx = 4, Ny = 4;
    dataType J1 = 1.0;
    dataType J2 = 0.0, dJ2 = 0.01;
    int J2_num = 0;
//     dataType J1 = {1.0,0.0};
//     dataType J2 = {0.0,0.0};
    switch (argc){
        case 3:
            Nx = atoi(argv[1]);
            Ny = atoi(argv[2]);
            break;
        case 4:
            Nx = atoi(argv[1]);
            Ny = atoi(argv[2]);
            J2 = dJ2 * atoi(argv[3]);
            break;
        case 5:
            Nx = atoi(argv[1]);
            Ny = atoi(argv[2]);
            J2_num = atoi(argv[3]);
            J2 = dJ2 * J2_num;
            nev = atoi(argv[4]);
            break;
    }
    int workerID, workerNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    #ifdef OMP_
        #pragma omp parallel
        {
            #pragma omp master
            ompThreadsNum = omp_get_num_threads();
        }
        if (workerID==MPI_MASTER) std::cout<<"openMP turned on with "<<ompThreadsNum<<" threads"<<std::endl;
    #else
        if (workerID==MPI_MASTER) std::cout<<"openMP turned off"<<std::endl;
    #endif
    
    for (int J2_num = 0; J2_num < 1; J2_num++){
        std::ofstream outfile;
        // data directory
        std::string dataDir = PROJECT_DATA_PATH+"/"+std::to_string(Nx)+"by"+std::to_string(Ny)+"/J2_"+std::to_string(J2_num);
        if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
        J2 = dJ2 * J2_num;
        if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;

    /*
        ****************************
        * Hamiltonian Construction *
        ****************************
    */
        TriAngLattice Lattice(Nx,Ny);
        int siteDim = 2;
        int dimList[] = {Nx*Ny/2, Nx*Ny/2};

        auto tic = std::chrono::system_clock::now();
        Basis B(siteDim, dimList);
        B.genTransReps(&Lattice);
//         idx_t count = 0;
//         for (idx_t j = 0; j < B.repCount.size(); j++) count += Lattice.N/B.repCount[j];
//         if(workerID==MPI_MASTER) std::cout<<"Total Rep number:"<<B.repIndList.size()<<". count = "<<count<<std::endl;
//         count = 0;
//         for (int i = 0; i < Lattice.N; ++i){
//             B.genKReps(i, &Lattice);
//             count += B.subDim;
//             if (workerID==MPI_MASTER) std::cout<<"subspace kInd = "<<i<<", size = "<<B.subDim<<std::endl;
//         }
        idx_t totDim = 0;
        for (int kIndex = 0;  kIndex < Lattice.N; kIndex++){
//         int kIndex = 2;
        B.genKReps(kIndex, &Lattice);
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"WorkerID:"<<workerID<<", subspace kInd = "<<kIndex<<", size = "<<B.subDim<<std::endl<<"********************"<<std::endl;
        totDim += B.subDim;
        HeisenbergOneHalf TriAngHeis(B.subDim, J1, J2);
        TriAngHeis.genSubMatMap(kIndex, &Lattice, &B, 2, 3);
//         if (workerID==MPI_MASTER) std::cout<<"Total size = "<<count<<std::endl;
        
//         HeisenbergOneHalf TriAngHeis(B.totDim, J1, J2);
//     //     TriAngHeis.genDiagMat(&Lattice, &B, 0, 3);
//         TriAngHeis.genMat(&Lattice, &B, 2, 3);
        auto toc = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = toc - tic;    
    //     TriAngHeis.saveAll(dataDir);

       std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<TriAngHeis.endRow - TriAngHeis.startRow<<"/"<<TriAngHeis.dim<<". Construction time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;
        
        MPI_Barrier(MPI_COMM_WORLD);

    /*
        *********************************
        * Diagonalization Using PARPACK *
        *********************************
    */   
        if (workerID==MPI_MASTER) std::cout<<"Initialize PARPACK..."<<std::endl;

//         PARPACKRealSolver<dataType> PDiag(&TriAngHeis, nev);
        // set ncv
        // PDiag.setNcv(5);
        PARPACKComplexSolver<double> PDiag(&TriAngHeis, nev);

        if (workerID==MPI_MASTER) std::cout<<"PARPACK Initialized."<<std::endl<<"Begin PARPACK Iteration and timer started..."<<std::endl;
        tic = std::chrono::system_clock::now();

        PDiag.run();

        toc = std::chrono::system_clock::now();
        elapsed_seconds = toc - tic;
        if (workerID==MPI_MASTER) std::cout<<"INFO:"<<PDiag.info_<<". Total iteration:"<<PDiag.iparam_[2]<<". Total time:"<<elapsed_seconds.count()*1000<<" milliseconds"<<std::endl;

        if (workerID==MPI_MASTER) std::cout<<"Begin post processing..."<<std::endl;

        PDiag.postRun();

        toc = std::chrono::system_clock::now();
        elapsed_seconds = toc - tic;
        if (workerID==MPI_MASTER) std::cout<<"INFO:"<<PDiag.info_<<". Total ev found:"<<PDiag.iparam_[4]<<". Total time:"<<elapsed_seconds.count()*1000<<" milliseconds"<<std::endl;
        if (workerID==MPI_MASTER) {
            std::cout<<"Eigenvalues: ";
            for (int i = 0; i < nev; ++i) std::cout<<PDiag.d_pt[i]<<", ";
            std::cout<<std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
//         }
        
//         // save eigen values
//         if (workerID==MPI_MASTER){
//             save<dataType>(PDiag.d_pt, nev, &outfile, dataDir + "/eigval");
//         }


    /*
        ******************************
        * Static Spin-Spin Structure *
        ******************************
    */
//         if (workerID==MPI_MASTER) std::cout<<"Begin SS ..."<<std::endl;
//         dataType* vecTmp1 = new dataType[PDiag.nlocmax_];
//         dataType* vecBuf1;
//         switch(PARTITION){
//             case ROW_PARTITION:
//                 vecBuf1 = new dataType[PDiag.ntot_];
//                 break;
//             case COL_PARTITION:
//                 vecBuf1 = new dataType[PDiag.nlocmax_];
//                 break;
//             default:break;
//         }

//         SSOneHalf SS(B.totDim); 
//         cdouble val;
//         std::vector<cdouble> ssvals;
//         for (int i = 0; i < Lattice.N; ++i){
//             SS.genPairMat(&Lattice, &B, i);
//             SS.MxV(PDiag.V_pt, vecTmp1, vecBuf1,PARTITION);
//             mpiDot<dataType, dataType>(PDiag.V_pt, vecTmp1, &val, PDiag.nloc_);
//             //if (workerID==MPI_MASTER) std::cout<<"SS"<<i<<" = "<<std::setprecision(8)<<val<<std::endl;
//             ssvals.push_back(val);
//             SS.clear();
//         }
//         // save ss(i)
//         if (workerID==MPI_MASTER) save<cdouble>(ssvals.data(), Lattice.N, &outfile, dataDir + "/spinspin");
//         delete [] vecBuf1; 
//         delete [] vecTmp1;
//         toc = std::chrono::system_clock::now();
//         elapsed_seconds = toc - tic;    
//         if (workerID==MPI_MASTER) std::cout<<"SS time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;

    /*
        *******************************
        * Dynamic Spin-Spin Structure *
        *******************************
    // */
//         if (workerID==MPI_MASTER) std::cout<<"Begin Sqw ..."<<std::endl;
//         cdouble* vecTmp = new cdouble[PDiag.nlocmax_];
//         cdouble* gstate = new cdouble[PDiag.nlocmax_];
//         cdouble* vecBuf;
//         switch(PARTITION){
//             case ROW_PARTITION:
//                 vecBuf = new cdouble[PDiag.ntot_];
//                 break;
//             case COL_PARTITION:
//                 vecBuf = new cdouble[PDiag.nlocmax_];
//                 break;
//             default:break;
//         }
//         if (workerID==MPI_MASTER) std::cout<<"Begin dynamic spin-spin structure calculation..."<<std::endl;
//         SzqOneHalf Szq(B.totDim);
//         if (workerID==MPI_MASTER) std::cout<<"ntot:"<<Szq.ntot<<". nlocmax:"<<Szq.nlocmax<<". nloc:"<<Szq.nloc<<std::endl;
//         for (int i = 1; i < Lattice.N; ++i){
//             BasisXY q; 
//             q[0] = Lattice.KLattice_[i].coordxy[0];
//             q[1] = Lattice.KLattice_[i].coordxy[1];
//             Szq.genMat(&Lattice, &B, q);
//     //         cdouble* vecTmp = new cdouble[Szq.nlocmax];
//     //         cdouble* gstate = new cdouble[Szq.nlocmax];
//             for (idx_t j = 0; j < Szq.nloc; j++) gstate[j] = (cdouble) PDiag.V_pt[j];
//         //     cdouble *vecBuf = new cdouble[Szq.nlocmax];
//             Szq.MxV(gstate, vecTmp, vecBuf, PARTITION);
//             cdouble val;
//             mpiDot(gstate, vecTmp, &val, Szq.nloc);
//     //         if (workerID==MPI_MASTER) std::cout<<"Ground state <Szq> = "<<val<<std::endl;

//             int krylovDim = 250;
//             LANCZOSIterator<dataType>lanczos(&TriAngHeis, krylovDim);
//             lanczos.run(vecTmp);
//     //         if (workerID==MPI_MASTER){
//     //             std::cout<<std::endl;
//     //             for (int i = 0; i < krylovDim; ++i) std::cout<<"a["<<i<<"] = "<<lanczos.alpha[i]<<std::endl;
//     //             std::cout<<std::endl;
//     //             for (int i = 0; i < krylovDim - 1; ++i) std::cout<<"b["<<i<<"] = "<<lanczos.beta[i]<<std::endl;
//     //             std::cout<<std::endl;
//     //         }

//         // save alpha, beta
//             if (workerID==MPI_MASTER){
//                 save<cdouble>(lanczos.alpha.data(), krylovDim, &outfile, dataDir + "/alpha_q" + std::to_string(i));
//                 save<cdouble>(lanczos.beta.data(), krylovDim-1, &outfile, dataDir + "/beta_q" + std::to_string(i));
//             }
//             Szq.clear();
//         }

//         // delete allocated memory
//         delete [] gstate;
//         delete [] vecBuf; 
//         delete [] vecTmp;
//         toc = std::chrono::system_clock::now();
//         elapsed_seconds = toc - tic;    
//         if (workerID==MPI_MASTER) std::cout<<"Sqw time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;
//         MPI_Barrier(MPI_COMM_WORLD);
        }
        if (workerID==MPI_MASTER) std::cout<<"Total Dim:"<<B.totDim<<". Sum of subspece dim:"<<totDim<<std::endl;
    }
    
/*
    *******
    * End *
    *******
*/
    MPI_Finalize();

    return 0;
}

