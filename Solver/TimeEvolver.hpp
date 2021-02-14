//
// TimeEvolver.hpp
// ED
//
// Created by tatang on 11/19/19
//
//

#ifndef TimeEvolver_hpp
#define TimeEvolver_hpp

#include "LANCZOSIterator.hpp"
#include "../Algebra/algebra.hpp"

#include <mpi.h>
#include <vector>
#include <complex>
#include <cmath>
#include "omp.h"

/*
    Krylov Subspace Algorithm for exp(-iH)*v
    Appropriate for sparse matrix H.
    usually krydim~10-30 should be accurate enough
*/
template <class T>
class TimeEvolver: public LANCZOSIterator<T>{
private:
    T* vec;
    std::vector<cdouble> tmp1, tmp2;
    std::vector<double> U;
public:
    TimeEvolver(){};
    TimeEvolver(T* initVec, BaseMatrix<dataType> *M, int krydim):LANCZOSIterator<T>(M,krydim,ALPHA_BETA_Q),vec(initVec),tmp1(krydim),tmp2(krydim),U(krydim*krydim){};
    ~TimeEvolver(){};

    T* getVec(){return vec;}
    // time evolve one step
    void evolve(double expFac);
};

// vec = exp(i*expFac*H)*vec
template <class T>
void TimeEvolver<T>::evolve(double expFac){
    LANCZOSIterator<T>::run(vec);
    // Tri = Z*Diag*Z^
    MKL::diagTri(&(LANCZOSIterator<T>::alpha), &(LANCZOSIterator<T>::beta), &U);
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::krylovDim; ++i){
        tmp1[i] = std::exp(-CPLX_I*expFac*LANCZOSIterator<T>::alpha.at(i))*U[i*LANCZOSIterator<T>::krylovDim];
        tmp2[i] = 0.0;
    } 
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::krylovDim; ++i){
        for (int j = 0; j < LANCZOSIterator<T>::krylovDim; j++){
            tmp2[i] += U[i+j*LANCZOSIterator<T>::krylovDim] * tmp1[j];
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < LANCZOSIterator<T>::M_->getnloc(); ++i){
        cdouble sum = 0.0;
        for(int j = 0; j < LANCZOSIterator<T>::krylovDim; j++){
            sum += LANCZOSIterator<T>::Q[i + j *LANCZOSIterator<T>::M_->getnloc()] * tmp2[j];
        }
        vec[i] = sum;
    }
}
#endif //TimeEvolver_hpp
