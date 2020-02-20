#include "TimeEvolver.hpp"

int main(int argc, const char * argv[]){
    MPI_Init(NULL, NULL);
    ind_int dim = 10000;
    int krydim =100;
    double expFac = 0.01;
    SparseMatrix<cdouble> sMat(dim);
    ind_int counter = 0;
    sMat.rowColInitList.push_back(counter);
    for (ind_int i = 0; i < dim; i++){
        sMat.valList.push_back(1.0*i);
        sMat.colList.push_back(i);
        counter++;
        sMat.rowColInitList.push_back(counter);
    }
    std::vector<cdouble> vec, vecf;
    for (ind_int i = 0; i < dim; i++){
        vec.push_back(1.0);
        vecf.push_back(std::exp(CPLX_I*expFac*i)*vec.at(i));
    }
    TimeEvolver<cdouble> test(vec.data(), &sMat, krydim);
    test.runStep(expFac);
    double err = 0.0;
    for (ind_int i = 0; i < krydim; i++){
        err += std::abs(vec.at(i)*test.vecNorm-vecf.at(i));
        // if (delta>err) err=delta;
    }
    std::cout<<"err bound:"<<err<<std::endl;
    return 0;
}
