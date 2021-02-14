#include <iostream>
#include "Utils/utils.hpp"
#include "Utils/paras.hpp"
#include "Global/config.hpp"
#include "Solver/PARPACKSolver.hpp"
#include "Solver/TimeEvolver.hpp"

using namespace std;

int main(){

    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);

    if (workerID == MPI_MASTER) {
        path.make();
        path.print();
        modelPara.print(path.parameterFile);
        pulsePara.print(path.pumpFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    setBasics(modelPara, latt, Bi, H);
    // cout<<setprecision(10);
    PARPACKComplexSolver<double> PDiag(H.get(), 1);
    PDiag.diag();
    cdouble* gstate = PDiag.getEigvec();

    // time evolution
    setpulse(pulsePara, pulse);
    H->setPeierls(&pulse);
    H->construct();
    
    int krylovDim = 15;
    TimeEvolver<cdouble> Tevol(gstate, H.get(), krylovDim);

    Nocc occ(latt.get(), Bi.get()); 
    occ.genMat();

    while (H->next()) {
        Tevol.evolve(-0.01);
        occ.count(Tevol.getVec());
    }
    if (workerID == MPI_MASTER) {
        occ.save(path.pumpDir);
    }

    MPI_Finalize();
    return 0;
}