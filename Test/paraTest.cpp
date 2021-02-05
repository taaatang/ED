#include <iostream>
#include "Utils/paras.hpp"

int main(){
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
    Parameters para("../Input/config.txt");
    std::ofstream outfile("para.txt");
    std::cout<<"Hello from "<<workerID<<".\n";
    MPI_Barrier(MPI_COMM_WORLD);
    if (workerID == MPI_MASTER) {
        std::cout<<"print by "<<workerID<<".\n";
        para.print(outfile);
    }
    MPI_Finalize();
    return 0;
}