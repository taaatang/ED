//
//  main.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "operator/operators.hpp"
#include "solver/parpackSolver.hpp"
#include "solver/spectra.hpp"
#include "utils/timer.hpp"

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
    Timer timer;
    a_int nev = 3;
    int N, Nx=4, Ny=1, Nu=3, Nd=2;
    int kIndex = 3; // 0 is Gamma Point
    int PGRepIndex = -1;
    int rowPerThread = 1;
    int rowCount = 50;
    int rowPerIt = 1000;
    std::string label;
    double tdp_val=1.0, tpp_val=0.5, tppz_val=0.3, Vd=0.0, Vp=2.8, Vpz=3.0, Ud=8.5, Up=4, Udp=0.6;
    
    // infile<std::string>({&label},"../Input/current_label.txt");
    // infile<int>({&Nx, &Ny, &Nu, &Nd}, "../Input/lattice_input.txt");
    // infile<int>({&kIndex, &PGRepIndex}, "../Input/symm_input.txt");
    // infile<double>({&tdp_val,&tpp_val,&tppz_val,&Vd,&Vp,&Vpz,&Ud,&Up,&Udp},"../Input/Hubbard_params.txt");
    // infile<int>({&rowCount, &rowPerIt}, "Input/genmatBuf_input.txt");
    N = Nx * Ny;

    // data directory
    // std::string subDir = "sqOcta_Udp_"+std::to_string(Nx) + "x" + std::to_string(Ny)+"_"+std::to_string(Nu)+"u"+std::to_string(Nd)+"d";
    // std::string subDir = std::to_string(N);
    // std::string subDir = "/sq"+std::to_string(N)+"_"+std::to_string(Nu)+"u"+std::to_string(Nd)+"d";
    // std::string dataDirP = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/Conductivity_"+label;
    std::string dataDir = "/Volumes/Shared/ProjectData/Photodoping/data/sqOcta/4x1/3u2d/k3/sigma";
    std::ofstream outfile;
    // basis data path 
    bool BASIS_IS_SAVED = false;
    // std::string basisDir = PROJECT_DATA_PATH+"/" + subDir + "/kSpace/basis/"+std::to_string(kIndex);
    // std::string basisfile = basisDir + "/basis";
    // std::string normfile = basisDir + "/norm";
 
/*
    ************************************
    * Lattice and Basis Initialization *
    ************************************
*/
    // geometry class
    SquareLattice Lattice(Nx,Ny);
    Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    Lattice.addOrb({ORBITAL::Py,5,{0.0,-0.5,0.0}});
    Lattice.construct();
    
    int siteDim = 2;
    VecI occList{Nu, Nd};
    idx_t fullDim=0, totDim=0;
    // std::cout<<"Nu:"<<Nu<<", Nd:"<<Nd<<"\n";
    
    /*
        **********************
        * Basis Construction *
        **********************
    */
    
    timer.tik();
    Basis B(LATTICE_MODEL::HUBBARD, &Lattice, occList, kIndex, PGRepIndex);
    if(workerID==MPI_MASTER)std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED){
        // #ifdef DISTRIBUTED_BASIS
        // B.construct(basisfile, normfile, workerID, workerNum);
        // #else
        // B.construct(basisfile, normfile);
        // #endif
    } else {
        B.construct();
    }
    timer.tok();
    if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
    fullDim = B.getTotDim(); totDim += B.getSubDim();
    if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kInd ="<<kIndex<<", size="<<B.getLocDim()<<"/"<<B.getSubDim()<<"/"<<B.getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl; 
    
    // int num = 2;
    // for(int idx_vpz=0; idx_vpz<num; idx_vpz++){
    //     for(int idx_tpz=0; idx_tpz<num; idx_tpz++){
            // std::string dataDir = dataDirP+"/vpz_"+std::to_string(idx_vpz)+"_tpz_"+std::to_string(idx_tpz);
            // std::string dataDir = dataDirP;
            // if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
            /*
                ****************************
                * Hamiltonian Construction *
                ****************************
            */
            // double tdp_val=1.13, tpp_val=0.49, tppz_val=0.3, Vd=0.0, Vp=3.24, Vpz=3.0, Ud=8.5, Up=4.1;
            // double tdp_val=1, tpp_val=0.5, tppz_val=0.0+idx_tpz*0.5/(num-1), Vd=0.0, Vp=3.2, Vpz=2.0+idx_vpz*(3.2-2.0)/(num-1), Ud=8.5, Up=4;
            Link<dataType> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -tdp_val); tdpx.addLinkVec({0.5,0.0,0.0});
            Link<dataType> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, tdp_val); tdpy.addLinkVec({0.0,0.5,0.0});
            Link<dataType> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, tdp_val); tpxd.addLinkVec({0.5,0.0,0.0});
            Link<dataType> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -tpp_val); tpxpy.addLinkVec({0.5,0.5,0.0}).addLinkVec({-0.5,-0.5,0.0});
            Link<dataType> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, tpp_val); tpxpyp.addLinkVec({0.5,-0.5,0.0}).addLinkVec({-0.5,0.5,0.0});
            Link<dataType> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -tdp_val); tpyd.addLinkVec({0.0,0.5,0.0});

            bool NotConst=true, isOrdered=false;
            Link<dataType> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpxpz.addLinkVec({-0.5,0.0,0.5});
            Link<dataType> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -tppz_val, NotConst, isOrdered); tpzpx.addLinkVec({0.5,0.0,0.5});
            Link<dataType> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpypz.addLinkVec({0.0,-0.5,0.5});
            Link<dataType> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -tppz_val, NotConst, isOrdered); tpzpy.addLinkVec({0.0,0.5,0.5});

            Link<dataType> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpxpzp.addLinkVec({0.5,0.0,0.5});
            Link<dataType> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, tppz_val, NotConst, isOrdered); tpzpxp.addLinkVec({-0.5,0.0,0.5});
            Link<dataType> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpypzp.addLinkVec({0.0,0.5,0.5});
            Link<dataType> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, tppz_val, NotConst, isOrdered); tpzpyp.addLinkVec({0.0,-0.5,0.5});

            Link<dataType> ndnpx(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Px}, Udp); ndnpx.addLinkVec({0.5,0.0,0.0}).addLinkVec({-0.5,0.0,0.0});
            Link<dataType> ndnpy(LINK_TYPE::HUBBARD_U, {ORBITAL::Dx2y2, ORBITAL::Py}, Udp); ndnpy.addLinkVec({0.0,0.5,0.0}).addLinkVec({0.0,-0.5,0.0});
            
            timer.tik();
            Hubbard<dataType> H(&Lattice, &B, 1);
            H.pushLinks({tdpx,tdpy,tpxd,tpxpy,tpxpyp,tpyd,tpxpz,tpxpzp,tpzpx,tpzpxp,tpypz,tpypzp,tpzpy,tpzpyp,ndnpx,ndnpy});
            // H.pushLinks({tpxpz,tpxpzp,tpzpx,tpzpxp,tpypz,tpypzp,tpzpy,tpzpyp});
            H.pushV({ORBITAL::Dx2y2},Vd).pushV({ORBITAL::Px,ORBITAL::Py},Vp).pushV({ORBITAL::Pzu, ORBITAL::Pzd},Vpz);
            H.pushU({ORBITAL::Dx2y2},Ud).pushU({ORBITAL::Px,ORBITAL::Py,ORBITAL::Pzu, ORBITAL::Pzd},Up);
            if(workerID==MPI_MASTER) std::cout<<"begin H gen..."<<std::endl;
            #ifdef DISTRIBUTED_BASIS
                H.construct(rowCount,rowPerIt);
            #else
                H.construct();
            #endif
            timer.tok();

            if(workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.getnloc()<<"/"<<H.getDim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
                <<". Construction time:"<<timer.elapse()<<"ms."<<std::endl;
        /*
            *********************************
            * Diagonalization Using PARPACK *
            *********************************
        */
            PARPACKComplexSolver<double> PDiag(&H, nev);
            PDiag.diag();
            cdouble* w0 = PDiag.getEigval();
            dataType* gstate = PDiag.getEigvec();
            MPI_Barrier(MPI_COMM_WORLD);

            // Nocc occ(&Lattice, &B); 
            // occ.genMat();
            // for(int state=0;state<nev;state++){
            //     dataType* vecpt = PDiag.getEigvec(state);
            //     double Nd = occ.count(ORBITAL::Dx2y2, vecpt);
            //     double Npx = occ.count(ORBITAL::Px, vecpt);
            //     double Npy = occ.count(ORBITAL::Py, vecpt);
            //     if(workerID==MPI_MASTER)std::cout<<"state "<<state<<": Nd "<<Nd<<", Npx "<<Npx<<", Npy "<<Npy<<"\n";
            // }


        /*
            ****************
            * Conductivity *
            * **************
        */
            std::vector<std::string> lables{"x", "y", "z"};
            for (const auto& label:lables) {
                timer.tik();
                Current J(&Lattice, &B, label);
                if(label=="z"){
                    J.pushLinks({tpxpz,tpxpzp,tpzpx,tpzpxp,tpypz,tpypzp,tpzpy,tpzpyp}); // Jz
                }else if(label=="y"){
                    J.pushLinks({tdpy,tpyd,tpxpy,tpxpyp,tpypz,tpypzp,tpzpy,tpzpyp}); // Jy
                }else if(label=="x"){
                    J.pushLinks({tdpx,tpxd,tpzpx,tpzpxp,tpxpz,tpxpzp,tpxpy,tpxpyp}); // Jx
                }else{
                    exit(1);
                }
                std::cout<<"begin construct J...\n";
                #ifdef DISTRIBUTED_BASIS
                    J.construct(rowCount, rowPerIt);
                #else
                    J.construct();
                #endif
                std::cout<<"init spectra\n";
                int krylovDim=400;
                SPECTRASolver<dataType> spectra(&H, w0[0], &J, gstate, H.getDim(), krylovDim);
                std::cout<<"compute spectra..\n";
                spectra.compute();
                // save alpha, beta
                std::cout<<"saving spectra...\n";
                if (workerID==MPI_MASTER){
                    std::string dataPath = dataDir + "/" + label;
                    spectra.saveData(dataPath);
                } 
                timer.tok();
                if (workerID==MPI_MASTER) std::cout<<"SigmaW time:"<<timer.elapse()<<" milliseconds."<<std::endl<<std::endl;
            }
    //     }
    // }

/*
    *******
    * End *
    *******
*/
    MPI_Finalize();
    return 0;
}