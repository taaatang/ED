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

    std::ofstream outfile("modelPara.txt");
    if (workerID == MPI_MASTER) {
        modelPara.print(outfile);
    }
    setBasics(modelPara, latt, Bi, H);
    // cout<<setprecision(10);
    PARPACKComplexSolver<double> PDiag(H.get(), 1);
    PDiag.diag();
    cdouble* gstate = PDiag.getEigvec();

    setpulse(pulsePara, pulse);
    std::string TimePath = PROJECT_DATA_PATH + "/sqOcta/4x1/2u2d/k2/pump/x/" + tostr(pulse.getw());
    mkdir_fs(TimePath);
    H->setPeierls(&pulse);
    H->construct();
    TimeEvolver<cdouble> Tevol(gstate, H.get(), 15);
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
        std::ofstream outfile;
        save<double>(dx2y2.data(),(int)dx2y2.size(),&outfile,TimePath+"/dx2y2");
        save<double>(px.data(),(int)px.size(),&outfile,TimePath+"/px");
        save<double>(py.data(),(int)py.size(),&outfile,TimePath+"/py");
        save<double>(pzu.data(),(int)pzu.size(),&outfile,TimePath+"/pzu");
        save<double>(pzd.data(),(int)pzd.size(),&outfile,TimePath+"/pzd");
    }
    MPI_Finalize();
    return 0;
}