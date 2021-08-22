#include <vector>
#include <iostream>
#include "algebra/algebra.hpp"

using namespace std;
using val_t = double;
using Vec = std::vector<val_t>;

template <typename T>
void print(string name, std::vector<T>& v, bool cond = false) {
    if (cond) cout<<name<<", size:"<<v.size()<<", cap:"<<v.capacity()<<"\nval:"<<v<<"\n";
}

int main ( ) {
    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpiInfo(workerID, workerNum);
    bool isMaster = (workerID == 0);
    Vec x{1.0, 2.0, 3.0};
    auto size = x.size();
    Vec y(size, 0.0);
    MKL::copy(x.size(), x.data(), y.data());
    print("y", y, isMaster);
    MKL::scale(size, 2.0, x.data());
    print("2*x", x, isMaster);
    MKL::axpy(size, 2.0, y.data(), x.data());
    print("y+x", x, isMaster);

    if (isMaster) {
        cout<<"x.dot(x):"<<MKL::dot(size, x.data(), x.data())<<"\n";
        cout<<"norm(x):"<<MKL::norm(size, x.data())<<"\n";
    }
    auto dotres = MKL::mpiDot(size, x.data(), x.data());
    auto normres = MKL::mpiNorm(size, x.data());
    if (isMaster) {
        cout<<"mpi x.dot(x):"<<dotres<<"\n";
        cout<<"mpi norm(x):"<<normres<<"\n";
    }
    MPI_Finalize();
    return 0;
}