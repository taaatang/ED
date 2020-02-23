//
//  main.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "globalPara.hpp"
#include "HelperClass.hpp" // Timer
#include "Geometry.hpp"
#include "Basis.hpp"
#include "Operators.hpp"
#include "utils.hpp"
#include "PARPACKSolver.hpp"
#include "LANCZOSIterator.hpp"
#include "TimeEvolver.hpp"
#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
#include <stdlib.h> // system
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
    Timer timer;
    a_int nev;
    int N1, N2, N;
    int kIndex = -1; // Gamma Point
    double spin = 0.5;
    int rowPerThread = 1;
    bool BASIS_IS_SAVED = false;
    std::ifstream infile("ground_input.txt");
    infile>>N1>>N2>>nev;
    N = N1 * N2;
    int kStart = 0, kEnd = N;
    int J2Start = 0, J2End = 101, J2Step = 1;
    infile.close();
    double J1 = 1.0;
    double J2 = 0.0, dJ2 = 0.01;
    // data directory
    std::string subDir = std::to_string(N1) + "by" + std::to_string(N2);
    // std::string subDir = std::to_string(N);
    std::string dataDirP = PROJECT_DATA_PATH+"/"+subDir+"/kSpace";
/*
    ********************
    * MPI and OMP Info *
    ********************
*/
    int workerID, workerNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &workerID);
    MPI_Comm_size(MPI_COMM_WORLD, &workerNum);
    if (workerID==MPI_MASTER) std::cout<<"Total MPI Workers:"<<workerNum<<std::endl; 
    OMP_Info(workerID);    
/*
    ************************************
    * Lattice and Basis Initialization *
    ************************************
*/
    // geometry class
    TriAngLattice Lattice(N1,N2);
    Lattice.addOrb({});
    // int N1 = 2, N2 = 2;
    // N = N1 * N2;
    // SquareLattice Lattice(N1,N2);
    // Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    // Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    Lattice.construct();
    // if (workerID==MPI_MASTER) Lattice.print();
    int siteDim = 2;
    VecI occList{N/2, N-N/2};
    ind_int fullDim=0, totDim=0;
    // scan kIndex
    for (int kIndex = kStart;  kIndex < kEnd; kIndex++){
        std::ofstream outfile;
        std::string basisDir = PROJECT_DATA_PATH+"/" + subDir + "/kSpace/Basis/"+std::to_string(kIndex);
        std::string basisfile = basisDir + "/basis";
        std::string normfile = basisDir + "/norm";
        /*
            **********************
            * Basis Construction *
            **********************
        */
        timer.tik();
        Basis B(LATTICE_MODEL::HEISENBERG, &Lattice, occList, kIndex);
        std::cout<<"begin construc basis..."<<std::endl;
        if (BASIS_IS_SAVED) B.gen(basisfile, normfile);
        else B.gen();
        timer.tok();
        fullDim = B.getTotDim(); totDim += B.getSubDim();
        if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kInd ="<<kIndex<<", size="<<B.getSubDim()<<"/"<<B.getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl;
        if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
        /*
            ****************************
            * Hamiltonian Construction *
            ****************************
        */
        /*
            **************
            * HEISENBERG *
            **************
        */
        Link<double> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0);
        Link<double> J2Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, false);
        J1Link.addLinkVec(VecD{1.0,0.0,0.0}).addLinkVec(VecD{0.0,1.0,0.0}).addLinkVec(VecD{1.0,-1.0,0.0});
        J2Link.addLinkVec(VecD{1.0,1.0,0.0}).addLinkVec(VecD{-1.0,2.0,0.0}).addLinkVec(VecD{2.0,-1.0,0.0});
        timer.tik();
        Heisenberg<dataType> H(&Lattice, &B, 2);
        H.pushLink(J1Link).pushLink(J2Link);
        H.genMatPara(rowPerThread);
        timer.tok();
        /*
            ***********
            * HUBBARD *
            ***********
        */
        // double tdp_val=1.13, tpp_val=0.49, tppz_val=0.3, Vd=0.0, Vp=3.24, Vpz=3.0, Ud=8.5, Up=4.1;
        // Link<double> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -tdp_val); tdpx.addLinkVec({0.5,0.0,0.0});
        // Link<double> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, tdp_val); tdpy.addLinkVec({0.0,0.5,0.0});
        // Link<double> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, tdp_val); tpxd.addLinkVec({0.5,0.0,0.0});
        // Link<double> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -tpp_val); tpxpy.addLinkVec({0.5,0.5,0.0}).addLinkVec({-0.5,-0.5,0.0});
        // Link<double> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, tpp_val); tpxpyp.addLinkVec({0.5,-0.5,0.0}).addLinkVec({-0.5,0.5,0.0});
        // Link<double> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -tdp_val); tpyd.addLinkVec({0.0,0.5,0.0});

        // Link<double> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, tppz_val, false); tpxpz.addLinkVec({-0.5,0.0,0.5});
        // Link<double> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -tppz_val, false); tpxpzp.addLinkVec({0.5,0.0,0.5});
        // Link<double> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -tppz_val, false); tpzpx.addLinkVec({0.5,0.0,0.5});
        // Link<double> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, tppz_val, false); tpzpxp.addLinkVec({-0.5,0.0,0.5});
        // Link<double> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, tppz_val, false); tpxpz.addLinkVec({0.0,-0.5,0.5});
        // Link<double> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -tppz_val, false); tpypzp.addLinkVec({0.0,0.5,0.5});
        // Link<double> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -tppz_val, false); tpzpx.addLinkVec({0.0,0.5,0.5});
        // Link<double> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, tppz_val, false); tpzpxp.addLinkVec({0.0,-0.5,0.5});

        // timer.tik();
        // if (workerID==MPI_MASTER){J1Link.print(); J2Link.print();}
        // Hubbard<dataType> H(&Lattice, &B, 1);
        // H.pushLink(tdpx).pushLink(tdpy).pushLink(tpxd).pushLink(tpxpy).pushLink(tpxpyp).pushLink(tpyd);
        // H.pushLink(tpxpz).pushLink(tpxpzp).pushLink(tpzpx).pushLink(tpzpxp).pushLink(tpypz).pushLink(tpypzp).pushLink(tpzpy).pushLink(tpzpyp);
        // tdpx.print();tdpy.print();tpxd.print();tpxpy.print();tpxpyp.print();tpyd.print();
        // H.pushV({ORBITAL::Dx2y2},Vd).pushV({ORBITAL::Px,ORBITAL::Py},Vp);
        // H.pushU({ORBITAL::Dx2y2},Ud).pushU({ORBITAL::Px,ORBITAL::Py},Up);
        // H.printV();H.printU();
        // H.genMatPara(rowPerThread);
        // timer.tok();

        if(workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.get_nloc()<<"/"<<H.get_dim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
            <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

        timer.tik();
        SSOp<dataType> SS(&Lattice,&B);
        SS.genMatPara();
        timer.tok();
        if(workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local SS dimension:"<<SS.get_nloc()<<"/"<<SS.get_dim()<<", Local Hamiltonian non-zero elements count:"<<SS.nzCount()\
            <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;


        // scan J2
        for (int J2_num = J2Start; J2_num < J2End; J2_num += J2Step){
            if (workerID==MPI_MASTER) std::cout<<"**********************"<<std::endl<<"Begin J2 = "<<J2<<std::endl<<"*************************"<<std::endl;
            J2 = dJ2 * J2_num;
            std::ofstream outfile;
            // data directory
            std::string dataDir = dataDirP+"/eigs/J2_"+std::to_string(J2_num);
            if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
            // set J2 parameter
            H.setVal(J2Link,J2);
            MPI_Barrier(MPI_COMM_WORLD);
        /*
            *********************************
            * Diagonalization Using PARPACK *
            *********************************
    */
            PARPACKComplexSolver<double> PDiag(&H, nev);

            // project into total spin subspace
            // std::vector<cdouble> initVec(H.get_nloc());
            // randInit<cdouble>(initVec);
            // SS.project(spin,initVec.data());
            // PDiag.setStartVec(initVec.data());
            // MPI_Barrier(MPI_COMM_WORLD);
            // PDiag.diag(spin, &SS);
            PDiag.diag();

            std::vector<cdouble> stot;
            for (int i = 0; i < nev; i++){
                stot[i] = SS.vMv(PDiag.getEigvec(i),PDiag.getEigvec(i));
                if(workerID==MPI_MASTER)std::cout<<i<<"th eigen state total spin s(s+1) = "<<stot<<std::endl;
            }
            // save eigen values
            if (workerID==MPI_MASTER){
                save<dataType>(PDiag.getEigval(), nev, &outfile, dataDir + "/eigval_k" + std::to_string(kIndex));
                save<cdouble>(stot.data(), nev, &outfile, dataDir + "/eigval_spin_k" + std::to_string(kIndex));
            }
            MPI_Barrier(MPI_COMM_WORLD);
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