//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
//

#ifndef LANCZOSIterator_hpp
#define LANCZOSIterator_hpp

#include "../Global/globalPara.hpp"
#include "../Operator/SparseMatrix.hpp"
#include "../Utils/utils.hpp"
#include <mpi.h>
#include <vector>
#include <complex>
#include <cmath>
#include "mkl.h"


template <class T>
class LANCZOSIterator{
protected:
    // local vectorsize
    BaseMatrix<T> *M_;
    LANCZOS_OPTION option;
    int krylovDim;
    // stores alpha, beta
    std::vector<double> alpha, beta;
    cdouble vecNorm;
    // q_j-1, q_j, r
    cdouble *q_pre, *q, *r;
    // AQ = QH + residual
    cdouble *Q; 
public:
    LANCZOSIterator(){};
    LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt=LANCZOS_DEFAULT_OPTION);
    ~LANCZOSIterator();
    
    void init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt);
    void run(cdouble* vec);
    // reorthogonalize
    void runReOth(cdouble* vec);
};

template <class T>
LANCZOSIterator<T>::LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
    init(M, krydim, opt);
}

template <class T>
void LANCZOSIterator<T>::init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
    M_ = M;
    option = opt;
    krylovDim = krydim;
    q_pre = new cdouble[M_->getnlocmax()];
    q = new cdouble[M_->getnlocmax()];
    r = new cdouble[M->getnlocmax()];
    if (option == ALPHA_BETA_Q){
        Q = new cdouble[krydim * M_->getnloc()];
    }
}

template <class T>
LANCZOSIterator<T>::~LANCZOSIterator(){
    delete [] q_pre;
    delete [] q;
    delete [] r;
    if (option == ALPHA_BETA_Q) delete [] Q;
}

template <class T>
void LANCZOSIterator<T>::run(cdouble* vec){
    alpha.clear();
    beta.clear();
    cdouble result;
    // q = vec/|vec|
    cdouble zero(0.0,0.0), one(1.0,0.0);
    veqv(q, vec, M_->getnloc());
    vConjDotv(vec, vec, &result, M_->getnloc());
    vecNorm = std::sqrt(result);
    if (std::norm(result) < INFINITESIMAL){
        std::cout<<"norm(vecIn) = "<<std::norm(result)<<". Check if the input vec is 0."<<std::endl;
        return;
    }
    vdeva(q, vecNorm, M_->getnloc());
    // Q1 = [q]
    if (option == ALPHA_BETA_Q) veqv(&Q[0], q, M_->getnloc());
    // r = Aq
    M_->MxV(q, r);
    // a1 = q^r
    vConjDotv(q, r, &result, M_->getnloc());
    alpha.push_back(std::real(result));
    // r = r - a1*q
    saxpy(r, -(cdouble)alpha[0], q, M_->getnloc());
    // beta[1] = |r|
    vConjDotv(r, r, &result, M_->getnloc());
    beta.push_back(std::real(std::sqrt(result)));
    if (M_->getWorkerID()==MPI_MASTER && beta[0]==0.0){
        std::cout<<"LANZOSIteration stops for beta[0] = 0."<<std::endl;
        return;
    }
    
    for (int i = 1; i < krylovDim; ++i){
        // q_pre = q
        veqv(q_pre, q, M_->getnloc());
        // q = r/beta[i-1]
        veqatv(q, one/beta[i-1], r, M_->getnloc());
        //Qi = [Qi-1,q]
        if (option == ALPHA_BETA_Q) veqv(&Q[i*M_->getnloc()], q, M_->getnloc());
        // r = Aq
        M_->MxV(q, r);
        //r = r - beta[i-1]*q_pre
        saxpy(r, -(cdouble)beta[i-1], q_pre, M_->getnloc());
        //alpha[i] = q^r
        vConjDotv(q, r, &result, M_->getnloc());
        alpha.push_back(std::real(result));
        // r = r - alpha[i]*q
        saxpy(r, -(cdouble)alpha[i], q, M_->getnloc());
        // beta[i] = |r|
        vConjDotv(r, r, &result, M_->getnloc());
        beta.push_back(std::real(std::sqrt(result)));
        if (M_->getWorkerID()==MPI_MASTER && beta[i]==zero){
            std::cout<<"LANZOSIteration stops for beta["<<i<<"] = 0."<<std::endl;
            return;
        }
    }
}

#endif
