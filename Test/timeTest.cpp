#include "../operator/operators.hpp"
#include "../solver/parpackSolver.hpp"
#include "../solver/timeEvolver.hpp"
#include "../pulse/pulse.hpp"
#include "../utils/timer.hpp"

#include <fstream>
#include <iomanip>
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
    int workerID, workerNum;
    mpi_info(workerID, workerNum);

    Timer timer;
    a_int nev = 1;
    int kIndex = 0, PGRepIndex=-1;
    int Nx=2, Ny=2, Nu=2, Nd=2, N;
    int rowPerThread = 1;
    int widx = 15;
    // infile<int>({&widx},"Input/freq_input.txt");
    double w = 0.1*widx;
    // std::ifstream infile("time_input.txt");
    // infile.close();
    // infile>>Nx>>Ny;
    N = Nx * Ny;
    // data directory
    std::string subDir = "sqOcta"+std::to_string(Nx) + "by" + std::to_string(Ny)+"_"+std::to_string(Nu)+"u"+std::to_string(Nd)+"d";
    std::string dataDirP = PROJECT_DATA_PATH+"/"+subDir+"/kSpace/TimeEvolve/k"+std::to_string(kIndex);

    // if (workerID==MPI_MASTER){
    //     system(("mkdir -p " + dataDirP).c_str());

    // } 

    bool BASIS_IS_SAVED = false;  
/*
    ************************************
    * Lattice and Basis Initialization *
    ************************************
*/
    // geometry class
    SquareLattice Lattice(Nx,Ny);
    Lattice.addOrb({ORBITAL::Dx2y2,0,{0.0,0.0,0.0}}).addOrb({ORBITAL::Px,1,{0.5,0.0,0.0}}).addOrb({ORBITAL::Py,2,{0.0,0.5,0.0}});
    Lattice.addOrb({ORBITAL::Pzu,3,{0.0,0.0,0.5}}).addOrb({ORBITAL::Pzd,4,{0.0,0.0,-0.5}});
    // Lattice.addOrb({ORBITAL::Py,5,{0.0,-0.5,0.0}});
    Lattice.construct();

    int siteDim = 2;
    VecI occList{Nu, Nd};
    idx_t fullDim=0, totDim=0;
    std::ofstream outfile;
    std::string basisDir = PROJECT_DATA_PATH+"/" + subDir + "/kSpace/basis/"+std::to_string(kIndex);
    std::string basisfile = basisDir + "/basis";
    std::string normfile = basisDir + "/norm";
    /*
        **********************
        * Basis Construction *
        **********************
    */
    timer.tik();
    Basis B(MODEL::HUBBARD, &Lattice, occList, kIndex, PGRepIndex);
    if(workerID==MPI_MASTER)std::cout<<"begin construc basis..."<<std::endl;
    if (BASIS_IS_SAVED){
        #ifdef DISTRIBUTED_BASIS
        B.construct(basisfile, normfile, workerID, workerNum);
        #else
        B.construct(basisfile, normfile);
        #endif
    } else {
        B.construct();
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
    std::cout<<std::setprecision(10);
    // double tdp_val=1.13, tpp_val=0.49, tppz_val=0.3, Vd=0.0, Vp=3.24, Vpz=3.0, Ud=8.5, Up=4.1;
    double tdp_val=1, tpp_val=0.5, tppz_val=0.3, Vd=0.0, Vp=2.8, Vpz=3.0, Ud=8.5, Up=4;
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
    if(workerID==MPI_MASTER) std::cout<<"begin H gen..."<<std::endl;
    #ifdef DISTRIBUTED_BASIS
        H.construct(rowCount,rowPerIt);
    #else
        H.construct();
    #endif
    timer.tok();
    if(workerID==MPI_MASTER) std::cout<<"WorkerID:"<<workerID<<". Local Hamiltonian dimension:"<<H.getnloc()<<"/"<<H.getDim()<<", Local Hamiltonian non-zero elements count:"<<H.nzCount()\
        <<". Construction time:"<<timer.elapse()<<" milliseconds."<<std::endl;

    // std::ofstream outfile;
    // // data directory
    // std::string dataDir = dataDirP+"/eigs/J2_"+std::to_string(J2_num);
    // if (workerID==MPI_MASTER) system(("mkdir -p " + dataDir).c_str());
    MPI_Barrier(MPI_COMM_WORLD);
    int fluNum = 2;
    double pow_start = -1.0, pow_end = 3.0;
    for (int fluidx = 0; fluidx < fluNum; fluidx++){
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
        // double Fluence = std::pow(10.0,pow_start+(pow_end-pow_start)/(fluNum-1)*fluidx);
        double Fluence = 1.5;
        std::string TimePath = dataDirP + "/w_"+std::to_string(widx)+"_flu_"+std::to_string(fluidx)+"_";

        int krydim = 15;
        double dt = 0.01;
        int tsteps = 15000;
        double a = 0.2429; //nm
        double width = 10.0;
        
        Pulse pulse(w,width,dt,tsteps);
        pulse.seta(a);
        pulse.setFluence(Fluence);
        pulse.print();
        TimeEvolver<cdouble> Tevol(gstate, &H, krydim);
        Nocc occ(&Lattice, &B); 
        occ.genMat();
        timer.tik();
        VecD dx2y2,px,py,pzu,pzd;
        dx2y2.push_back(occ.count(ORBITAL::Dx2y2, Tevol.getVec()));
        px.push_back(occ.count(ORBITAL::Px,Tevol.getVec()));
        py.push_back(occ.count(ORBITAL::Py,Tevol.getVec()));
        pzu.push_back(occ.count(ORBITAL::Pzu,Tevol.getVec()));
        pzd.push_back(occ.count(ORBITAL::Pzd,Tevol.getVec()));

        for(int tstep = 0; tstep < tsteps; tstep++){
            double Aa = pulse.getAa(); pulse.next();
            cdouble factor = std::exp(CPLX_I*Aa);
            H.setVal(1,factor);
            H.setVal(2,std::conj(factor));
            Tevol.evolve(-dt);

            dx2y2.push_back(occ.count(ORBITAL::Dx2y2, Tevol.getVec()));
            px.push_back(occ.count(ORBITAL::Px,Tevol.getVec()));
            py.push_back(occ.count(ORBITAL::Py,Tevol.getVec()));
            pzu.push_back(occ.count(ORBITAL::Pzu,Tevol.getVec()));
            pzd.push_back(occ.count(ORBITAL::Pzd,Tevol.getVec()));

            if(tstep%150==0){
                if(workerID==MPI_MASTER) std::cout<<"\ntime step: "<<tstep<<"\n"\
                    <<"d:"<<occ.count(ORBITAL::Dx2y2,Tevol.getVec())<<"\n"\
                    <<"px:"<<occ.count(ORBITAL::Px,Tevol.getVec())<<"\n"\
                    <<"py:"<<occ.count(ORBITAL::Py,Tevol.getVec())<<"\n"\
                    <<"pzu:"<<occ.count(ORBITAL::Pzu,Tevol.getVec())<<"\n"\
                    <<"pzd:"<<occ.count(ORBITAL::Pzd,Tevol.getVec())<<"\n";
            }
        }
        // std::ofstream outfile;
        // save<double>(dx2y2.data(),(int)dx2y2.size(),&outfile,TimePath+"dx2y2");
        // save<double>(px.data(),(int)px.size(),&outfile,TimePath+"px");
        // save<double>(py.data(),(int)py.size(),&outfile,TimePath+"py");
        // save<double>(pzu.data(),(int)pzu.size(),&outfile,TimePath+"pzu");
        // save<double>(pzd.data(),(int)pzd.size(),&outfile,TimePath+"pzd");
        timer.tok();
        if(workerID==MPI_MASTER) std::cout<<"total iterations:"<<tsteps<<". time:"<<timer.elapse()<<"ms."<<std::endl;
    }
    MPI_Finalize();
    return 0;
}
