#include <iostream>
#include "utils/paras.hpp"
#include "global/config.hpp"
#include "solver/parpackSolver.hpp"
#include "solver/timeEvolver.hpp"

using namespace std;

int main(){
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
    std::ofstream outfile("para.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    if (workerID == MPI_MASTER) {
        para.print(outfile);
    }
    setlatt(para, latt);
    setbasis(para, Bi, latt.get());
    Bi->construct();
    setham(para, H, latt.get(), Bi.get());
    H->printLinks();
    H->construct();
    // if (workerID == MPI_MASTER) {
    //     latt->print();
    //     Bi->print();
    // }
    cout<<setprecision(10);
    PARPACKComplexSolver<double> PDiag(H.get(), 1);
    PDiag.diag();
    cdouble* gstate = PDiag.getEigvec();

    setpulse(pulsePara, pulse);
    H->setPeierls(&pulse);
    H->construct();
    TimeEvolver<cdouble> Tevol(gstate, H.get(), 15);
    Nocc occ(latt.get(), Bi.get()); 
    occ.genMat();
    int tstep = 0;
    while (H->next()) {
        Tevol.evolve(-0.01);
        if(tstep%150==0){
            if(workerID==MPI_MASTER) std::cout<<"\ntime step: "<<tstep<<"\n"\
                        <<"d:"<<occ.count(ORBITAL::Dx2y2,Tevol.getVec())<<"\n"\
                        <<"px:"<<occ.count(ORBITAL::Px,Tevol.getVec())<<"\n"\
                        <<"py:"<<occ.count(ORBITAL::Py,Tevol.getVec())<<"\n"\
                        <<"pzu:"<<occ.count(ORBITAL::Pzu,Tevol.getVec())<<"\n"\
                        <<"pzd:"<<occ.count(ORBITAL::Pzd,Tevol.getVec())<<"\n";
        }
        ++tstep;
    }
    MPI_Finalize();
    return 0;
}