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

    bool COMPUTE_SS = false;
    bool COMPUTE_SQW = false;
    bool COMPUTE_RAMAN = true;
    a_int nev = 5;
    // int Nx = 6, Ny = 6;
    // int N = Nx * Ny;
    LATTICE_MODEL model = LATTICE_MODEL::HEISENBERG;
    int Nx, Ny, N=12, Nu=6, Nd=6;
    int kIndex = 0; // Gamma Point
    int PGRepIndex = -1;
    double spin = 0.5;
    double t1 = 1.0, t2 = 0.2, J1= 1.0, J2 = 0.1, Jk = 0.0;
    double dJ = 0.01;
    int J2idx, Jkidx;

    Timer timer;

    int rowPerThread = 1;

    bool BASIS_IS_SAVED = false;

    // infile<int>({&N, &Nu, &Nd}, "../Input/lattice_input.txt");
    // infile<int>({&kIndex, &PGRepIndex}, "../Input/symm_input.txt");
    // infile<int>({&J2idx, &Jkidx}, "../Input/params_input.txt");
    for(J2idx=0;J2idx<51;J2idx++){
        for(Jkidx=0;Jkidx<51;Jkidx++){
    J2 = dJ * double(J2idx);
    Jk = dJ * double(Jkidx);

    // data directory
    // std::string subDir = std::to_string(Nx) + "by" + std::to_string(Ny);
    // std::string subDir = "N"+tostr(N)+"Nu"+tostr(Nu)+"Nd"+tostr(Nd);
    std::string subDir = tostr(N)+"_test";
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
    Basis *B = new Basis(model, &Lattice, occList, kIndex);
    if (workerID==MPI_MASTER)std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED) B->gen(basisfile, normfile);
    else B->gen();
    timer.tok();
    if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kIdx ="<<kIndex<<", PGidx = "<<PGRepIndex<<", size="<<B->getSubDim()<<"/"<<B->getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl;
    if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
    
    // VecD J2s{0.0, 0.16, 0.3, 0.6};
    // VecD Jks{0.0, 0.2, 0.4, 0.8};
    // for(auto J2 : J2s){
    //     for(auto Jk : Jks){
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
    Heisenberg<dataType>* H = new Heisenberg<dataType>(&Lattice, B, 1);
    H->pushLinks({J1Link, J2Link, JkLink});
    // HtJ<dataType> H(&Lattice, B, 1);
    // H.pushLinks({t1Link, t2Link, J1Link, J2Link});
    H->genMatPara();
    delete B;
    timer.tok();
    std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H->get_nloc()<<"/"<<H->get_dim()<<", Local Hamiltonian non-zero elements count:"<<H->nzCount()\
            <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

// for(int J2_num = 0; J2_num<101; J2_num++){
//     J2 = dJ2 * J2_num;
    std::ofstream outfile;
    // data directory
    std::string dataDir = PROJECT_DATA_PATH+"/"+ subDir +"/kSpace/Spectra/J2_"+tostr(J2idx)+"_Jk_"+tostr(Jkidx);
    if (workerID==MPI_MASTER) mkdir_fs(dataDir);
    if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;
    // H.setVal(1,J2);
    /*
        *********************************
        * Diagonalization Using PARPACK *
        *********************************
    */
    PARPACKComplexSolver<double> PDiag(H, nev);
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
    double wgs = ws[0].real();
    for (int state_idx = 0; state_idx<nev; state_idx++){
        if((ws[state_idx]).real()<wgs) wgs = ws[state_idx].real();
    }
    for (int state_idx = 0; state_idx<nev; state_idx++){
        if(std::abs((ws[state_idx]).real()-wgs)<tol) gs_idx.push_back(state_idx);
    }
    if(workerID==MPI_MASTER) std::cout<<"Ground State deg:"<<gs_idx.size()<<"\n";
/*
    *************************************
    * Static Spin-Spin Structure Factor *
    *************************************
*/
    if (COMPUTE_SS){
        Basis *B = new Basis(model, &Lattice, occList, kIndex);
        if (BASIS_IS_SAVED) B->gen(basisfile, normfile);
        else B->gen();
        timer.tik();
        SSOp<dataType> SS(&Lattice,B);
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
        delete B;
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
        ind_int Hdim = H->get_dim();
        H->clearBuf();
        delete H;
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Sqw ..."<<std::endl<<"********************"<<std::endl;
        for (int kIndexf = 0; kIndexf < Lattice.getSiteNum(); kIndexf++){
            timer.tik();
            if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<kIndexf<<std::endl<<"********************"<<std::endl;
            std::string basisDirp = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/Basis/"+std::to_string(kIndexf);
            std::string basisfilep = basisDirp + "/basis";
            std::string normfilep = basisDirp + "/norm";
            Basis *B = new Basis(model, &Lattice, occList, kIndex);
            if (BASIS_IS_SAVED) B->gen(basisfile, normfile);
            else B->gen();
            Basis* Bp = new Basis(model, &Lattice, occList, kIndexf);
            if (BASIS_IS_SAVED) {Bp->gen(basisfilep, normfilep);
            }else{Bp->gen();}
            // <Bp|Szq|B>, q = k_B - k_Bp
            SzkOp<dataType> Szq(&Lattice, B, Bp);
            Szq.genMatPara();

            delete B;

            Heisenberg<dataType> Hp(&Lattice, Bp, 1);
            Hp.pushLinks({J1Link,J2Link,JkLink});
            Hp.genMatPara();

            delete Bp;

            // Hp.setVal(1, J2);
            timer.tok();
            if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<Hp.get_nloc()<<"/"<<Hp.get_dim()<<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
            for(auto idx : gs_idx){
                cdouble w0 = ws[idx];
                dataType* state = PDiag.getEigvec(idx);
                SPECTRASolver<dataType> spectra(&Hp, w0, &Szq, state, Hdim, krylovDim);
                spectra.compute();
                // save alpha, beta
                if (workerID==MPI_MASTER){
                    std::string dataPath = dataDir + "/sqw_k" + std::to_string(kIndex) + "_kp" + std::to_string(kIndexf)+"_state_"+tostr(idx);
                    spectra.saveData(dataPath);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            } 

            Hp.clearBuf();

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

        // // A1
        // Link<dataType> A1Link_1(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1/2.0);
        // A1Link_1.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
        // Link<dataType> A1Link_2(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J2/2.0*3.0);
        // A1Link_2.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
        // std::vector<Link<dataType>> A1Links{A1Link_1, A1Link_2};
        // // E2_1
        // Link<dataType> E21Link_1(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1/2.0);
        // E21Link_1.addLinkVec({1.0,0.0,0.0});
        // Link<dataType> E21Link_2(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -J1/4.0);
        // E21Link_2.addLinkVec({0.0,1.0,0.0}).addLinkVec({-1.0,1.0,0.0});
        // Link<dataType> E21Link_3(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J2/4.0*3.0);
        // E21Link_3.addLinkVec({2.0,-1,0.0}).addLinkVec({1.0,1.0,0.0});
        // Link<dataType> E21Link_4(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -J2/2.0*3);
        // E21Link_4.addLinkVec({-1.0,2.0,0.0});
        // std::vector<Link<dataType>> E21Links{E21Link_1,E21Link_2,E21Link_3,E21Link_4};
        // // E2_2
        // Link<dataType> E22Link_1(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1/4.0*std::sqrt(3.0));
        // E22Link_1.addLinkVec({0.0,1.0,0.0});
        // Link<dataType> E22Link_2(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -J1/4.0*std::sqrt(3.0));
        // E22Link_2.addLinkVec({-1.0,1.0,0.0});
        // Link<dataType> E22Link_3(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J2/4.0*3.0*std::sqrt(3.0));
        // E22Link_1.addLinkVec({1.0,1.0,0.0});
        // Link<dataType> E22Link_4(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, -J2/4.0*3.0*std::sqrt(3.0));
        // E22Link_2.addLinkVec({2.0,-1.0,0.0});
        // std::vector<Link<dataType>> E22Links{E22Link_1,E22Link_2,E22Link_3,E22Link_4};
        // // A2
        // std::vector<VecD> A2vecs{{1.0,1.0,0.0},{1.0,0.0,0.0}, {1.0,0.0,0.0},{2.0,-1.0,0.0},\
        // {0.0,1.0,0.0},{1.0,0.0,0.0}, {1.0,0.0,0.0},{1.0,-1.0,0.0}, {1.0,1.0,0.0},{2.0,0.0,0.0}, {2.0,0.0,0.0},{2.0,-1.0,0.0}};
        // int vecNum = A2vecs.size();
        // for (int rot = 0; rot < 5; rot++){
        //     for (int vecidx = 0; vecidx < vecNum; vecidx++){
        //         int cur_idx = rot*vecNum+vecidx;
        //         VecD vec = Lattice.rotate(A2vecs.at(cur_idx));
        //         A2vecs.push_back(vec);
        //     }
        // }
        // Link<dataType> A2Link(LINK_TYPE::CHIRAL_K, {ORBITAL::SINGLE, ORBITAL::SINGLE, ORBITAL::SINGLE}, CPLX_I*4*J1*std::sqrt(3*J1*J2), true);
        // for(auto vec : A2vecs) A2Link.addLinkVec(vec);
        // std::vector<Link<dataType>> A2Links{A2Link};

        // std::vector<std::vector<Link<dataType>>> LinksList{A1Links,E21Links,E22Links,A2Links};
        // std::vector<std::string> RamanLabels{"A1","E21","E22","A2"};
        
        // for(int opidx=0;opidx<4;opidx++){
        //     H->clearBuf();
        //     Basis *B = new Basis(model, &Lattice, occList, kIndex);
        //     if (BASIS_IS_SAVED) B->gen(basisfile, normfile);
        //     else B->gen();
        //     Heisenberg<dataType> R(&Lattice, B);
        //     R.pushLinks(LinksList.at(opidx));
        //     R.genMatPara();
        //     delete B;

        //     timer.tok();
        //     if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Raman Operator construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
        //     for(auto idx : gs_idx){
        //         cdouble w0 = ws[idx];
        //         dataType* state = PDiag.getEigvec(idx);
        //         SPECTRASolver<dataType> spectra(H, w0, &R, state, H->get_dim(), krylovDim);
        //         spectra.compute();
        //         // save alpha, beta
        //         if (workerID==MPI_MASTER){
        //             std::string dataPath = dataDir + "/raman_k" + std::to_string(kIndex)+"_"+RamanLabels.at(opidx)+"_state_"+tostr(idx);
        //             spectra.saveData(dataPath);
        //         } 
        //         MPI_Barrier(MPI_COMM_WORLD);
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        // }

        VecD plzX{1.0,0.0,0.0}, plzY{-1.0,2.0,0.0};
        std::vector<VecD> plz{plzX,plzY};
        std::vector<std::string> plzLabel{"x","y"};
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin Raman ..."<<std::endl<<"********************"<<std::endl;
        timer.tik();
        if (workerID==MPI_MASTER) std::cout<<"********************"<<std::endl<<"Begin kIndex = "<<kIndex<<std::endl<<"********************"<<std::endl;
        for(int i = 0; i < 1; i++){
            for(int j = 0; j < plz.size(); j++){
                H->clearBuf();
                Basis *B = new Basis(model, &Lattice, occList, kIndex);
                if (BASIS_IS_SAVED) B->gen(basisfile, normfile);
                else B->gen();

                RamanOp<dataType> R(&Lattice, B);
                R.pushLinks({J1Link,J2Link});
                R.setplz(plz[i],plz[j]);
                R.genMatPara();

                delete B;
                
                timer.tok();
                if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Raman Operator construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;
                for(auto idx : gs_idx){
                    cdouble w0 = ws[idx];
                    dataType* state = PDiag.getEigvec(idx);
                    SPECTRASolver<dataType> spectra(H, w0, &R, state, H->get_dim(), krylovDim);
                    spectra.compute();
                    // save alpha, beta
                    if (workerID==MPI_MASTER){
                        std::string dataPath = dataDir + "/raman_k" + std::to_string(kIndex)+"_"+plzLabel[i]+plzLabel[j]+"_state_"+tostr(idx);
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
        // Heisenberg<dataType> Rc(&Lattice, B, 1);
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
    //     }
    // }
    
    if(!COMPUTE_SQW){
        delete H;
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

