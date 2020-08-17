//
//  SpectraMain.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//


#include "../Operator/Operators.hpp"
#include "../Solver/PARPACKSolver.hpp"
#include "../Solver/Spectra.hpp"
#include "../Utils/timer.hpp"

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
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
    mpi_info(workerID, workerNum);

    bool COMPUTE_SS = true;
    bool COMPUTE_SQW = true;
    bool COMPUTE_RAMAN = true;
    a_int nev = 1;
    // int Nx = 6, Ny = 6;
    // int N = Nx * Ny;
    LATTICE_MODEL model = LATTICE_MODEL::HEISENBERG;
    int Nx, Ny, N, Nu, Nd;
    int kIndex = -1; // Gamma Point
    int PGRepIndex = -1;
    double spin = 0.5;
    double t1 = 1.0, t2 = 0.2, J1= 1.0, J2 = 0.0, Jk = 0.0;
    double dJ2 = 0.01;

    Timer timer;

    int rowPerThread = 1;

    bool BASIS_IS_SAVED = false;

    infile<int>({&N, &Nu, &Nd}, "../Input/lattice_input.txt");
    infile<int>({&kIndex, &PGRepIndex}, "../Input/symm_input.txt");
    infile<double>({&J1, &J2, &Jk}, "../Input/params_input.txt");

    // data directory
    // std::string subDir = std::to_string(Nx) + "by" + std::to_string(Ny);
    // std::string subDir = "N"+tostr(N)+"Nu"+tostr(Nu)+"Nd"+tostr(Nd);
    std::string subDir = "N"+tostr(N)+"_test";
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
    VecI occList{Nu, Nd};
    Basis B(model, &Lattice, occList, kIndex);
    if (workerID==MPI_MASTER)std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED) B.gen(basisfile, normfile);
    else B.gen();
    timer.tok();
    if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kIdx ="<<kIndex<<", PGidx = "<<PGRepIndex<<", size="<<B.getSubDim()<<"/"<<B.getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl;
    if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
    
    VecD J2s{0.0, 0.16, 0.3, 0.6};
    VecD Jks{0.0, 0.2, 0.4, 0.8};
    for(auto J2 : J2s){
        for(auto Jk : Jks){
    /*
        ****************************
        * Hamiltonian Construction *
        ****************************
    */
    // Link<dataType> t1Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t1);
    // Link<dataType> t2Link(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t2);
    Link<dataType> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1);
    Link<dataType> J2Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J2, true);
    Link<dataType> JkLink(LINK_TYPE::CHIRAL_K, {ORBITAL::SINGLE, ORBITAL::SINGLE, ORBITAL::SINGLE},Jk, true);
    // t1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
    // t2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
    J1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
    J2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
    JkLink.addLinkVec({0.0,1.0,0.0}).addLinkVec({1.0,0.0,0.0}).addLinkVec({1.0,0.0,0.0}).addLinkVec({1.0,-1.0,0.0});
    timer.tik();
    Heisenberg<dataType> H(&Lattice, &B, 1);
    H.pushLinks({J1Link, J2Link, JkLink});
    // HtJ<dataType> H(&Lattice, &B, 1);
    // H.pushLinks({t1Link, t2Link, J1Link, J2Link});
    H.genMatPara();
    timer.tok();
    std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.get_nloc()<<"/"<<H.get_dim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
            <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

// for(int J2_num = 0; J2_num<101; J2_num++){
//     J2 = dJ2 * J2_num;
    std::ofstream outfile;
    // data directory
    std::string dataDir = PROJECT_DATA_PATH+"/"+ subDir +"/kSpace/Spectra/J2_"+tostr(J2)+"_Jk_"+tostr(Jk);
    if (workerID==MPI_MASTER) mkdir_fs(dataDir);
    if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;
    // H.setVal(1,J2);
    /*
        *********************************
        * Diagonalization Using PARPACK *
        *********************************
    */
    PARPACKComplexSolver<double> PDiag(&H, nev);
    PDiag.diag();
    // save eigen values
    if (workerID==MPI_MASTER){
        save<dataType>(PDiag.getEigval(), nev, &outfile, dataDir + "/eigval_k" + tostr(kIndex));
    }
    // eigval
    cdouble* ws = PDiag.getEigval();
    // dataType* state = PDiag.getEigvec();
    MPI_Barrier(MPI_COMM_WORLD);
    VecI gs_idx;
    double tol = 1e-8;
    for (int state_idx = 0; state_idx<nev; state_idx++){
        if(std::abs(ws[state_idx].real-ws[0].real)<tol) gs_idx.push_back(state_idx);
    }
    if(workerID==MPI_MASTER) std::cout<<"Ground State deg:"<<gs_idx.size()<<"\n";
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
            for(auto idx : gs_idx){
                dataType* state = PDiag.getEigvec(idx);
                SS.MxV(state, vecTmp.data());
                vConjDotv<dataType, dataType>(state, vecTmp.data(), &val, SS.get_nloc());
                val /= Lattice.getSiteNum();
                ssvals.push_back(val);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            
            if (workerID==MPI_MASTER) std::cout<<"SS "<<i<<" finished:"<<val<<std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        // save ss(i)
        if (workerID==MPI_MASTER) save<cdouble>(ssvals.data(), Lattice.getSiteNum(), &outfile, dataDir + "/spinspin_k"+std::to_string(kIndex));
        timer.tok();
        if (workerID==MPI_MASTER) std::cout<<"SS time:"<<timer.elapse()<<" milliseconds."<<std::endl;
    }
    

/*
    **************************************
    * Dynamic Spin-Spin Structure Factor *
    **************************************
*/
    if (COMPUTE_SQW){
        int krylovDim = 400;
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Sqw ..."<<std::endl<<"********************"<<std::endl;
        for (int kIndexf = 0; kIndexf < Lattice.getSiteNum(); kIndexf++){
            timer.tik();
            if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<kIndexf<<std::endl<<"********************"<<std::endl;
            std::string basisDirp = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/Basis/"+std::to_string(kIndexf);
            std::string basisfilep = basisDirp + "/basis";
            std::string normfilep = basisDirp + "/norm";
            Basis Bp = Basis(model, &Lattice, occList, kIndexf);
            if (BASIS_IS_SAVED) {Bp.gen(basisfilep, normfilep);
            }else{Bp.gen();}
            // <Bp|Szq|B>, q = k_B - k_Bp
            SzkOp<dataType> Szq(&Lattice, &B, &Bp);
            Szq.genMatPara();
            Heisenberg<dataType> Hp(&Lattice, &Bp, 1);
            Hp.pushLinks({J1Link,J2Link,JkLink});
            Hp.genMatPara();
            // Hp.setVal(1, J2);
            timer.tok();
            if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<Hp.get_nloc()<<"/"<<Hp.get_dim()<<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
            for(auto idx : gs_idx){
                cdouble w0 = ws[idx];
                dataType* state = PDiag.getEigvec(idx);
                SPECTRASolver<dataType> spectra(&Hp, w0, &Szq, state, H.get_dim(), krylovDim);
                spectra.compute();
                // save alpha, beta
                if (workerID==MPI_MASTER){
                    std::string dataPath = dataDir + "/sqw_k" + std::to_string(kIndex) + "_kp" + std::to_string(kIndexf)+"state_"+tostr(idx);
                    spectra.saveData(dataPath);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            } 
            timer.tok();
            if (workerID==MPI_MASTER) std::cout<<"Sqw time:"<<timer.elapse()<<" milliseconds."<<std::endl<<std::endl;
        }
    }
/*
    *********
    * Raman *
    *********
*/
    if (COMPUTE_RAMAN){
        int krylovDim = 400;
        VecD plzX{1.0,0.0,0.0}, plzY{1.0,std::sqrt(3),0.0};
        std::vector<VecD> plz{plzX,plzY};
        std::vector<std::string> plzLabel{"x","y"};
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Raman ..."<<std::endl<<"********************"<<std::endl;
        timer.tik();
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<kIndex<<std::endl<<"********************"<<std::endl;
        for(int i = 0; i < 1; i++){
            for(int j = 0; j < plz.size(); j++){
                RamanOp<dataType> R(&Lattice, &B);
                R.pushLinks({J1Link,J2Link});
                R.setplz(plz[i],plz[j]);
                R.genMatPara();
                timer.tok();
                if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Raman Operator construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
                for(auto idx : gs_idx){
                    cdouble w0 = ws[idx];
                    dataType* state = PDiag.getEigvec(idx);
                    SPECTRASolver<dataType> spectra(&H, w0, &R, state, H.get_dim(), krylovDim);
                    spectra.compute();
                    // save alpha, beta
                    if (workerID==MPI_MASTER){
                        std::string dataPath = dataDir + "/raman_k" + std::to_string(kIndex)+"_"+plzLabel[i]+plzLabel[j]+"state_"+tostr(idx);
                        spectra.saveData(dataPath);
                    } 
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }       
        timer.tok();
        if (workerID==MPI_MASTER) std::cout<<"Raman time:"<<timer.elapse()<<" milliseconds."<<std::endl<<std::endl;
        
        // timer.tik();
        // Heisenberg<dataType> Rc(&Lattice, &B, 1);
        // JkLink.setVal(1.0);
        // Rc.pushLinks({JkLink});
        // Rc.genMatPara();
        // timer.tok();
        // if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Chiral Raman Operator construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
        // MPI_Barrier(MPI_COMM_WORLD);
        // SPECTRASolver<dataType> spectra(&H, w0[0], &Rc, gstate, H.get_dim(), krylovDim);
        // spectra.compute();
        // // save alpha, beta
        // if (workerID==MPI_MASTER){
        //     std::string dataPath = dataDir + "/raman_k" + std::to_string(kIndex)+"_chiral";
        //     spectra.saveData(dataPath);
        // }
        // timer.tok();
        // if (workerID==MPI_MASTER) std::cout<<"Chiral Raman time:"<<timer.elapse()<<" milliseconds."<<std::endl<<std::endl;
    }
        }
    }
    

// }   
/*
    *******
    * End *
    *******
*/
    MPI_Finalize();
    return 0;
}

