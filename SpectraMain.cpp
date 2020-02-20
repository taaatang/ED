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
    bool COMPUTE_SS = true;
    bool COMPUTE_SQW = true;
    MPI_Init(NULL, NULL);
    int ompThreadsNum;
    a_int nev = 1;
    int Nx = 6, Ny = 6;
    int N = 36;
    dataType J1 = 1.0;
    dataType J2 = 0.0, dJ2 = 0.01;
    // int J2_num = 0;
    int kIndex = 0; // Gamma Point
    std::ifstream infile("input.txt");
    infile>>N>>kIndex;
    infile.close();
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
            // J2_num = atoi(argv[3]);
            break;
        case 5:
            Nx = atoi(argv[1]);
            Ny = atoi(argv[2]);
            // J2_num = atoi(argv[3]);
            kIndex = atoi(argv[4]);
            break;
    }
    // J2 = dJ2 * J2_num;
    // int N = Nx * Ny;
    
    for(int J2_num = 0; J2_num<101; J2_num++){
    J2 = dJ2 * J2_num;
    int workerID, workerNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl;
    
    OMP_Info(workerID);
    
    std::ofstream outfile;
    // data directory
    // std::string subDir = std::to_string(Nx) + "by" + std::to_string(Ny);
    std::string subDir = std::to_string(N);
    std::string dataDir = PROJECT_DATA_PATH+"/"+ subDir +"/kSpace/Spectra/J2_"+std::to_string(J2_num);
    if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
    
    if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;
    
    bool BASIS_IS_SAVED = false;
    std::string basisDir = PROJECT_DATA_PATH+"/" + subDir + "/kSpace/Basis/"+std::to_string(kIndex);
    std::string basisfile = basisDir + "/basis";
    std::string normfile = basisDir + "/norm";

    
/*
    ************************************
    * Lattice and Basis Initialization *
    ************************************
*/
    // TriAngLattice Lattice(Nx, Ny);
    TriAngLattice Lattice(N);
    int siteDim = 2;
    int dimList[] = {N/2, N-N/2};
    auto tic = std::chrono::system_clock::now();
    Basis B(MODEL, siteDim, dimList, &Lattice);
    if (BASIS_IS_SAVED) {
        B.gen(kIndex, basisfile, normfile);
    }else{
        B.gen(kIndex);
    }
    auto toc = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = toc - tic;
    
//         int kIndex = 2;
    if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"WorkerID:"<<workerID<<", subspace kInd = "<<kIndex<<", size = "<<B.subDim<<". Basis construct time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl<<"********************"<<std::endl;
//    std::cout<<"workerID:"<<workerID<<", B.kRep size:"<<B.indexList.size()<<", first and last elements:"<<B.indexList.at(0)<<" "<<B.indexList.at(B.subDim-1)<<std::endl;
//    std::cout<<"workerID:"<<workerID<<", B.norm size:"<<B.normList.size()<<", first and last elements:"<<B.normList.at(0)<<" "<<B.normList.at(B.subDim-1)<<std::endl;
/*
    ****************************
    * Hamiltonian Construction *
    ****************************
*/
    tic = std::chrono::system_clock::now();
    Heisenberg *TriAngHeis = new(std::nothrow) Heisenberg(&Lattice, &B, J1, J2); assert(TriAngHeis!=NULL);
    TriAngHeis->genSubMatMap(kIndex, 2, 3);
    toc = std::chrono::system_clock::now();
    elapsed_seconds = toc - tic;
//     TriAngHeis.saveAll(dataDir);

   if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<TriAngHeis->endRow - TriAngHeis->startRow<<"/"<<TriAngHeis->dim<<". Construction time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    ind_int ntot = TriAngHeis->get_ntot();
    ind_int nlocmax = TriAngHeis->get_nlocmax();

/*
    *********************************
    * Diagonalization Using PARPACK *
    *********************************
*/
    if (workerID==MPI_MASTER) std::cout<<"Initialize PARPACK..."<<std::endl;

//         PARPACKRealSolver<dataType> PDiag(&TriAngHeis, nev);
    // set ncv
    // PDiag.setNcv(5);
    PARPACKComplexSolver<double> *PDiag = new(std::nothrow) PARPACKComplexSolver<double>(TriAngHeis, nev); assert(PDiag!=NULL);

    if (workerID==MPI_MASTER) std::cout<<"PARPACK Initialized."<<std::endl<<"Begin PARPACK Iteration and timer started..."<<std::endl;
    tic = std::chrono::system_clock::now();

    PDiag->run();

    toc = std::chrono::system_clock::now();
    elapsed_seconds = toc - tic;
    if (workerID==MPI_MASTER) std::cout<<"INFO:"<<PDiag->info_<<". Total iteration:"<<PDiag->iparam_[2]<<". Total time:"<<elapsed_seconds.count()*1000<<" milliseconds"<<std::endl;

    if (workerID==MPI_MASTER) std::cout<<"Begin post processing..."<<std::endl;

    PDiag->postRun();

    toc = std::chrono::system_clock::now();
    elapsed_seconds = toc - tic;
    if (workerID==MPI_MASTER) std::cout<<"INFO:"<<PDiag->info_<<". Total ev found:"<<PDiag->iparam_[4]<<". Total time:"<<elapsed_seconds.count()*1000<<" milliseconds"<<std::endl;
    if (workerID==MPI_MASTER) {
        std::cout<<"Eigenvalues: ";
        for (int i = 0; i < nev; i++) std::cout<<PDiag->d_pt[i]<<", ";
        std::cout<<std::endl;
    }
    // Ground State
    dataType* gstate = new(std::nothrow) dataType[PDiag->nlocmax_]; assert(gstate!=NULL);
    double tol = 1e-8;
    int deg_count = 0;
    for (int eig_ind = 0; eig_ind < nev; eig_ind++){
        if (std::abs(PDiag->d_pt[eig_ind] - PDiag->d_pt[0]) > tol) continue;
        deg_count++;
        for (ind_int j = 0; j < PDiag->nloc_; j++) gstate[j] = (cdouble) PDiag->V_pt[eig_ind * PDiag->nloc_ + j];
    //    std::cout<<workerID<<"B subDim:"<<B.subDim<<" ground state nloc"<<PDiag->nloc_<<std::endl;
        
        delete TriAngHeis;
        delete PDiag;
        
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
        if (COMPUTE_SS){
            tic = std::chrono::system_clock::now();
            SSOneHalf *SS = new(std::nothrow) SSOneHalf(B.subDim); assert(SS!=NULL);
            if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin SS ..."<<std::endl<<"********************"<<std::endl;
            dataType* vecTmp1 = new(std::nothrow) dataType[SS->get_nlocmax()]; assert(vecTmp1!=NULL);
            dataType* vecBuf1;
            switch(PARTITION){
                case ROW_PARTITION:
                    vecBuf1 = new(std::nothrow) dataType[SS->get_ntot()]; assert(vecBuf1!=NULL);
                    break;
                case COL_PARTITION:
                    vecBuf1 = new(std::nothrow) dataType[SS->get_nlocmax()]; assert(vecBuf1!=NULL);
                    break;
                default:break;
            }
            cdouble val;
            vConjDotv(gstate, gstate, &val, SS->get_nloc());
            std::vector<cdouble> ssvals;
            for (int i = 0; i < Lattice.N; i++){
                val = 0.0;
        //        if (workerID==MPI_MASTER) std::cout<<"begin "<<i<<std::endl;
                SS->genPairMat(kIndex, &Lattice, &B, i);
        //        if (workerID==MPI_MASTER) std::cout<<"SS "<<i<<" created"<<std::endl;
                SS->MxV(gstate, vecTmp1, vecBuf1,PARTITION);
                vConjDotv<dataType, dataType>(gstate, vecTmp1, &val, SS->get_nloc());
                //if (workerID==MPI_MASTER) std::cout<<"SS"<<i<<" = "<<std::setprecision(8)<<val<<std::endl;
                val /= Lattice.N;
                ssvals.push_back(val);
                SS->clear();
                if (workerID==MPI_MASTER) std::cout<<"SS "<<i<<" finished:"<<val<<std::endl;
            }
            // save ss(i)
            if (workerID==MPI_MASTER) std::cout<<workerID<<" saving data"<<std::endl;
            if (workerID==MPI_MASTER) save<cdouble>(ssvals.data(), Lattice.N, &outfile, dataDir + "/spinspin_" + std::to_string(deg_count));
            delete [] vecBuf1;
            delete [] vecTmp1;
            delete SS;
            toc = std::chrono::system_clock::now();
            elapsed_seconds = toc - tic;
            if (workerID==MPI_MASTER) std::cout<<"SS time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        

    /*
        *******************************
        * Dynamic Spin-Spin Structure *
        *******************************
    */
    #define PRINT_
        if (COMPUTE_SQW){
            if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Sqw ..."<<std::endl<<"********************"<<std::endl;
    //        int numJ2 = 13
    //        int J2_List[numJ2] = {0, 4, 8, 12, 16, 18, 20, 24, 28, 32, 40, 60 80}
            assert(PARTITION==ROW_PARTITION);
            for (int i = 0; i < Lattice.N; i++){
                tic = std::chrono::system_clock::now();
                if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<i<<std::endl<<"********************"<<std::endl;
                std::string basisDirp = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/Basis/"+std::to_string(i);
                std::string basisfilep = basisDirp + "/basis";
                std::string normfilep = basisDirp + "/norm";
                Basis *Bp = new Basis(MODEL, siteDim, dimList, &Lattice);
                if (BASIS_IS_SAVED) {
                    Bp->gen(i, basisfilep, normfilep);
                }else{
                    Bp->gen(i);
                }
                
                SzqOneHalf Szq(Bp->subDim);
                // <Bp|Szq|B>, q = k_B - k_Bp
                BasisXY q;
                q[0] = Lattice.KLattice_[kIndex].coordxy[0] - Lattice.KLattice_[i].coordxy[0];
                q[1] = Lattice.KLattice_[kIndex].coordxy[1] - Lattice.KLattice_[i].coordxy[1];
                Szq.genMat(&Lattice, &B, Bp, q);
                MPI_Barrier(MPI_COMM_WORLD);

                cdouble* vecTmp = new cdouble[Szq.get_nlocmax()];
                cdouble* vecBuf;
                switch(PARTITION){
                    case ROW_PARTITION:
                        vecBuf = new(std::nothrow) cdouble[ntot]; assert(vecBuf!=NULL);
                        break;
                    case COL_PARTITION:
                        vecBuf = new(std::nothrow) cdouble[nlocmax]; assert(vecBuf!=NULL);
                        break;
                    default:break;
                }
                Szq.MxV(gstate, vecTmp, vecBuf, PARTITION);
                delete [] vecBuf;
        //        cdouble val;
        //        vConjDotv(gstate, vecTmp, &val, Szq.nloc);
        //        if (workerID==MPI_MASTER) std::cout<<"Ground state <Szq> = "<<val<<std::endl;
                Heisenberg *Hp = new Heisenberg(&Lattice, Bp, J1, J2);
                Hp->genSubMatMap(i, 2, 3);
                toc = std::chrono::system_clock::now();
                elapsed_seconds = toc - tic;
                if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<Hp->endRow - Hp->startRow<<"/"<<Hp->dim<<". Construction time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl;
                
                delete Bp;
                MPI_Barrier(MPI_COMM_WORLD);

                int krylovDim = 250;
                if (workerID==MPI_MASTER) std::cout<<"workerID:"<<workerID<<"begin lanczos ... "<<std::endl;
                LANCZOSIterator<dataType> lanczos(Hp, krylovDim);
                lanczos.run(vecTmp);
                // save alpha, beta
                if (workerID==MPI_MASTER){
                    save<cdouble>(&(lanczos.vecNorm), 1, &outfile, dataDir + "/vecNorm_" + std::to_string(deg_count) + "_k" + std::to_string(kIndex) + "_kp" + std::to_string(i));
                    save<double>(lanczos.alpha.data(), (int)lanczos.alpha.size(), &outfile, dataDir + "/alpha_" + std::to_string(deg_count) + "_k" + std::to_string(kIndex) + "_kp" + std::to_string(i));
                    save<double>(lanczos.beta.data(), (int)lanczos.beta.size(), &outfile, dataDir + "/beta_" + std::to_string(deg_count) + "_k" + std::to_string(kIndex) + "_kp" + std::to_string(i));
                }
                delete Hp;
                delete [] vecTmp;
                toc = std::chrono::system_clock::now();
                elapsed_seconds = toc - tic;
                if (workerID==MPI_MASTER) std::cout<<"Sqw time:"<<elapsed_seconds.count()*1000<<" milliseconds."<<std::endl<<std::endl;
            }
        }
    }
    // delete allocated memory
    delete [] gstate;
    }
   

    
/*
    *******
    * End *
    *******
*/
    MPI_Finalize();

    return 0;
}

