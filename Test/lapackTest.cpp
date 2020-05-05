#include <iostream>
#include <vector>
#include <cmath>
#include "mkl.h"
#include "omp.h"

int main(int argc, const char * argv[]) {
    int n = 10;
    std::vector<double> alpha;
    std::vector<double> beta;
    std::vector<double> Z;
    beta.resize(n-1,0.0);
    Z.resize(n*n);
    for (int i = 0; i < n; i++){
        alpha.push_back(i);
    }
    LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', n, alpha.data(), beta.data(), Z.data(), n);
    for (int i = 0; i < n; i++){
        std::cout<<"Eigen Value: "<<alpha.at(i)<<std::endl;
        std::cout<<"Eigen Vec: ";
        for (int j = i*n; j<(i+1)*n; j++) std::cout<<Z.at(j)<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    
}