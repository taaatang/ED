//
// TimeEvolver.hpp
// ED
//
// Created by tatang on 11/19/19
//
//

#ifndef TimeEvolver_hpp
#define TimeEvolver_hpp

#include "globalPara.hpp"
#include "SparseMatrix.hpp"
#include "LANCZOSIterator.hpp"
#include "utils.hpp"
#include "algebra.hpp"
#include <mpi.h>
#include <vector>
#include <complex>
#include <cmath>
//#include <Eigen/Eigenvalues>
#include "mkl.h"
#include "omp.h"

/*
    Krylov Subspace Algorithm for exp(-iH)*v
    Appropriate for sparse matrix H.
    usually krydim~10-30 should be accurate enough
*/
template <class T>
class TimeEvolver: public LANCZOSIterator<T>{
public:
    cdouble* vec;
    std::vector<cdouble> tmp1, tmp2;
    std::vector<double> U;
    TimeEvolver(){};
    TimeEvolver(cdouble* initVec, BaseMatrix<dataType> *M, int krydim):LANCZOSIterator<T>(M, krydim, ALPHA_BETA_Q){
        vec = initVec;
        tmp1.resize(krydim);
        tmp2.resize(krydim);
        U.resize(krydim*krydim);
    };
    ~TimeEvolver(){};

    // time evolve one step
    void runStep(double expFac);
};

template <class T>
void TimeEvolver<T>::runStep(double expFac){
    LANCZOSIterator<T>::run(vec);
    // Tri = Z*Diag*Z^
    diagTri(&(LANCZOSIterator<T>::alpha), &(LANCZOSIterator<T>::beta), &U);
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::krylovDim; i++){
        tmp1[i] = std::exp(CPLX_I*expFac*LANCZOSIterator<T>::alpha.at(i))*U[i*LANCZOSIterator<T>::krylovDim];
    } 
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::krylovDim; i++){
        for (int j = 0; j < LANCZOSIterator<T>::krylovDim; j++){
            tmp2[i] += U[i+j*LANCZOSIterator<T>::krylovDim] * tmp1[j];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::M_->get_nloc(); i++){
        cdouble sum = 0.0;
        for(int j = 0; j < LANCZOSIterator<T>::krylovDim; j++){
            sum += LANCZOSIterator<T>::Q[i + j *LANCZOSIterator<T>::M_->get_nloc()] * tmp2[j];
        }
        vec[i] = sum;
    }
}
#endif
