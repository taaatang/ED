#include <iostream>
#include "Utils/paras.hpp"
#include "Global/config.hpp"
#include "Solver/PARPACKSolver.hpp"

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
    Bi->gen();
    setham(para, H, latt.get(), Bi.get());
    H->printLinks();
    H->genMatPara();

    // if (workerID == MPI_MASTER) {
    //     latt->print();
    //     Bi->print();
    // }
    cout<<setprecision(10);
    PARPACKComplexSolver<double> PDiag(H.get(), 1);
    PDiag.diag();
    MPI_Finalize();
    return 0;
}