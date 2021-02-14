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
    VecD dx2y2,px,py,pzu,pzd;
    while (H->next()) {
        Tevol.evolve(-0.01);
        dx2y2.push_back(occ.count(ORBITAL::Dx2y2, Tevol.getVec()));
        px.push_back(occ.count(ORBITAL::Px,Tevol.getVec()));
        py.push_back(occ.count(ORBITAL::Py,Tevol.getVec()));
        pzu.push_back(occ.count(ORBITAL::Pzu,Tevol.getVec()));
        pzd.push_back(occ.count(ORBITAL::Pzd,Tevol.getVec()));
    }

    if (workerID == MPI_MASTER) {
        std::string pumpDir = path.pumpDir;
        save<double>(dx2y2.data(), (int)dx2y2.size(), pumpDir + "/dx2y2");
        save<double>(px.data(), (int)px.size(), pumpDir + "/px");
        save<double>(py.data(), (int)py.size(), pumpDir + "/py");
        save<double>(pzu.data(), (int)pzu.size(), pumpDir + "/pzu");
        save<double>(pzd.data(), (int)pzd.size(), pumpDir + "/pzd");
    }

    MPI_Finalize();
    return 0;
}