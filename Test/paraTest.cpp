#include <iostream>
#include "Utils/paras.hpp"
#include "Global/config.hpp"

int main(){
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
    std::ofstream outfile("para.txt");
    std::cout<<"Hello from "<<workerID<<".\n";
    MPI_Barrier(MPI_COMM_WORLD);
    if (workerID == MPI_MASTER) {
        std::cout<<"print by "<<workerID<<".\n";
        para.print(outfile);
    }
    setlatt(para, latt);
    setbasis(para, Bi, latt.get());
    if (workerID == MPI_MASTER) {
        latt->print();
    }
    MPI_Finalize();
    return 0;
}