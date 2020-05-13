#include "../Operator/Operators.hpp"
#include "../Solver/PARPACKSolver.hpp"
#include "../Solver/TimeEvolver.hpp"
#include "../Utils/timer.hpp"

#include <fstream>
#include <stdio.h>
#include <stdlib.h> // system
#include <mpi.h>
#include <omp.h>

int main(int argc, const char * argv[]){
    /*
    ****************************
    * Input And Initialization *
    ****************************
*/
    MPI_Init(NULL, NULL);
    Timer timer;
    a_int nev = 1;
    int kIndex = -1;
    int Nx=4, Ny=1, Nu=2, Nd=2, N;
    int rowPerThread = 1;
    // std::ifstream infile("time_input.txt");
    // infile.close();
    // infile>>N1>>N2;
    N = Nx * Ny;
    // data directory
    std::string subDir = "sqOcta"+std::to_string(Nx) + "by" + std::to_string(Ny);
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
    SquareLattice Lattice(N1,N2);
    Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    Lattice.addOrb({ORBITAL::Py,5,{0.0,-0.5,0.0}});
    Lattice.construct();

    int siteDim = 2;
    VecI occList{Nu, Nd};
    ind_int fullDim=0, totDim=0;
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
    Basis B(LATTICE_MODEL::HUBBARD, &Lattice, occList, kIndex, PGRepIndex);
    if(workerID==MPI_MASTER)std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED){
        #ifdef DISTRIBUTED_BASIS
        B.gen(basisfile, normfile, workerID, workerNum);
        #else
        B.gen(basisfile, normfile);
        #endif
    } else {
        B.gen();
    }
    timer.tok();
    if (workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". k-subspace Basis constructed:"<<timer.elapse()<<" milliseconds."<<std::endl;
    fullDim = B.getTotDim(); totDim += B.getSubDim();
    if (workerID==MPI_MASTER) std::cout<<std::endl<<"**********************"<<std::endl<<"Begin subspace kInd ="<<kIndex<<", size="<<B.getLocDim()<<"/"<<B.getSubDim()<<"/"<<B.getTotDim()<<std::endl<<"*************************"<<std::endl<<std::endl; 
    /*
        ****************************
        * Hamiltonian Construction *
        ****************************
    */
    /*
        ***********
        * HUBBARD *
        ***********
    */
    // double tdp_val=1.13, tpp_val=0.49, tppz_val=0.3, Vd=0.0, Vp=3.24, Vpz=3.0, Ud=8.5, Up=4.1;
    double tdp_val=1, tpp_val=0.5, tppz_val=0.5, Vd=0.0, Vp=3.2, Vpz=3.2, Ud=8.5, Up=4;
    Link<dataType> tdpx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px}, -tdp_val); tdpx.addLinkVec({0.5,0.0,0.0});
    Link<dataType> tdpy(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Py}, tdp_val); tdpy.addLinkVec({0.0,0.5,0.0});
    Link<dataType> tpxd(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Dx2y2}, tdp_val); tpxd.addLinkVec({0.5,0.0,0.0});
    Link<dataType> tpxpy(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, -tpp_val); tpxpy.addLinkVec({0.5,0.5,0.0}).addLinkVec({-0.5,-0.5,0.0});
    Link<dataType> tpxpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Py}, tpp_val); tpxpyp.addLinkVec({0.5,-0.5,0.0}).addLinkVec({-0.5,0.5,0.0});
    Link<dataType> tpyd(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Dx2y2}, -tdp_val); tpyd.addLinkVec({0.0,0.5,0.0});

    bool NotConst=false, isOrdered=true;
    Link<dataType> tpxpz(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpxpz.addLinkVec({-0.5,0.0,0.5});
    Link<dataType> tpxpzp(LINK_TYPE::HOPPING_T, {ORBITAL::Px, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpxpzp.addLinkVec({0.5,0.0,0.5});
    Link<dataType> tpzpx(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, -tppz_val, NotConst, isOrdered); tpzpx.addLinkVec({0.5,0.0,0.5});
    Link<dataType> tpzpxp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Px}, tppz_val, NotConst, isOrdered); tpzpxp.addLinkVec({-0.5,0.0,0.5});
    Link<dataType> tpypz(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, tppz_val, NotConst, isOrdered); tpypz.addLinkVec({0.0,-0.5,0.5});
    Link<dataType> tpypzp(LINK_TYPE::HOPPING_T, {ORBITAL::Py, ORBITAL::Pzu}, -tppz_val, NotConst, isOrdered); tpypzp.addLinkVec({0.0,0.5,0.5});
    Link<dataType> tpzpy(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, -tppz_val, NotConst, isOrdered); tpzpy.addLinkVec({0.0,0.5,0.5});
    Link<dataType> tpzpyp(LINK_TYPE::HOPPING_T, {ORBITAL::Pzd, ORBITAL::Py}, tppz_val, NotConst, isOrdered); tpzpyp.addLinkVec({0.0,-0.5,0.5});

    timer.tik();
    Hubbard<dataType> H(&Lattice, &B, 3);
    H.pushLinks({tdpx,tdpy,tpxd,tpxpy,tpxpyp,tpyd});
    H.pushLinks({tpxpz,tpxpzp,tpzpx,tpzpxp,tpypz,tpypzp,tpzpy,tpzpyp});
    H.pushV({ORBITAL::Dx2y2},Vd).pushV({ORBITAL::Px,ORBITAL::Py},Vp).pushV({ORBITAL::Pzu, ORBITAL::Pzd},Vpz);
    H.pushU({ORBITAL::Dx2y2},Ud).pushU({ORBITAL::Px,ORBITAL::Py,ORBITAL::Pzu, ORBITAL::Pzd},Up);
    H.genMatPara(&B);
    // H.genMat();
    timer.tok();
    if(workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.get_nloc()<<"/"<<H.get_dim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
        <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

    // std::ofstream outfile;
    // // data directory
    // std::string dataDir = dataDirP+"/eigs/J2_"+std::to_string(J2_num);
    // if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
    MPI_Barrier(MPI_COMM_WORLD);
/*
    *********************************
    * Diagonalization Using PARPACK *
    *********************************
*/
    PARPACKComplexSolver<double> PDiag(&H, nev);
    PDiag.diag();
    cdouble* gstate = PDiag.getEigvec();
/*
    ******************
    * Time Evolution *
    ******************
*/
    int krydim = 15;
    double dt = 0.01;
    int tsteps = 100;
    double t, t0=dt*(tsteps/2);
    double amp = 0.5, sigma = dt*(tsteps/5), freq = 1.5;
    TimeEvolver<cdouble> Tevol(gstate, &H, krydim);
    Nocc occ(&Lattice, &B); 
    occ.genMat();
    timer.tik();
    for(int tstep = 0; tstep < tsteps; tstep++){
        t = dt * tstep - t0;
        double A = GaussPulse(t,amp,sigma,freq);
        cdouble factor = std::exp(CPLX_I*A);
        H.setVal(1,factor);
        H.setVal(2,std::conj(factor));
        Tevol.evolve(-dt);
        if(workerID==MPI_MASTER) std::cout<<"time step: "<<t<<"\n"\
            <<"d:"<<occ.count(ORBITAL::Dx2y2,Tevol.getVec())<<"\n"\
            <<"px:"<<occ.count(ORBITAL::Px,Tevol.getVec())<<"\n"\
            <<"py:"<<occ.count(ORBITAL::Py,Tevol.getVec())<<"\n"\
            <<"pzu:"<<occ.count(ORBITAL::Pzu,Tevol.getVec())<<"\n"\
            <<"pzd:"<<occ.count(ORBITAL::Pzd,Tevol.getVec());
    }
    timer.tok();
    if(workerID==MPI_MASTER) std::cout<<"total iterations:"<<tsteps<<". time:"<<timer.elapse()<<"ms."<<std::endl;
    MPI_Finalize();
    return 0;
}
