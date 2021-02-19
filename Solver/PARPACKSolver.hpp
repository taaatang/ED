//
//  PARPACKSolver.hpp
//  ED
//
//  Created by tatang on 11/18/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#ifndef PARPACKSolver_hpp
#define PARPACKSolver_hpp
#include <assert.h>
#include <type_traits>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include "parpack.hpp"

#include "Global/globalPara.hpp"
#include "Utils/utils.hpp"
#include "Utils/timer.hpp"
#include "Operator/SparseMatrix.hpp"
#include "Operator/Operators.hpp"

template <class T>
class PARPACKSolver {
public:
    PARPACKSolver( ) { }
    PARPACKSolver(BaseMatrix<T> *M, a_int nev = 1);
    ~PARPACKSolver( ) { }

    T* getEigval( ) {return d_pt.data();}
    T* getEigvec(int ithVec=0) {return V_pt.data()+nloc_*ithVec;}

    void setStartVec(T* vec);
    void setNev(a_int nev);
    void setNcv(a_int ncv);
    void setMaxIter(int iterNum);
    void setRvec(a_int rvec);

    void diag( );
    void diag(double spin, SSOp<T>* SS);

private:
    void check( );

    void run( );
    void postRun( );

    // project to spin sector
    void run(double spin, SSOp<T>* SS);
    // project out states
    void run(T* states, int statesNum, double penalty=1000.0);

private:
    MPI_Fint MCW;
    BaseMatrix<T>* M_{nullptr};
    a_uint ntot_{0}, nlocmax_{0}, nloc_{0}, ldz_{0}, ldv_{0};
    int rank_;
    a_int ido_, info_, nev_{0}, ncv_{0}, lworkl_{0}, rvec_{0};
    double tol_{0.0};
    T sigma_;

    // output. ipntr(0) pointer to current operand vector in workd. 
    // ipntr(1) pointer to current result vector in workd
    std::array<a_int, 14> ipntr_;
    std::array<a_int, 11> iparam_;
    std::vector<T> workd_pt, workl_pt, V_pt, d_pt, z_pt, resid_pt;
    std::vector<T> rwork_pt, workev_pt;
    std::vector<a_int> select_pt;
};

template <class T>
PARPACKSolver<T>::PARPACKSolver(BaseMatrix<T> *H, a_int nev){
    if (H->getDim() <= 0) {
        std::cout<<"Empty Matrix for PARPACK Solver!\n";
        exit(1);
    }
    MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    // Input:use random initial vector. Output:0-normal exit,1-maximum iteration taken
    info_ = 0;
    // reverse communication flag  
    ido_ = 0; 

    // pointer to the hamitonian matrix
    M_ = H; 
    ntot_ = M_->getntot();
    // maximum local vector size
    nlocmax_ = M_->getnlocmax();
    // actual local vector size
    nloc_ = M_->getnloc();
    // number of eigenvalues to be computed
    nev_ = nev; 
    // number of Lanczos vectors generated at each iteration. number of colunm of V
    ncv_ = (2 * nev_ + 1) > PARPACK_MINNCV?(2 * nev_ + 1):PARPACK_MINNCV;
    assert_msg(nev_ > 0, "nev_ <= 0");
    assert_msg(ncv_ > 0, "ncv_ <= 0");
    if (ncv_ > (a_int)M_->getDim()) ncv_ = (a_int)M_->getDim()-1;
    ldz_ = nlocmax_ + 1;
    // leading dimension of V
    ldv_ = nlocmax_;
    rvec_ = 1;

    if constexpr (std::is_same<T, double>::value) {
        lworkl_ = ncv_ * (ncv_ + 8);
    } else if constexpr (std::is_same<T, cdouble>::value) {
        // dimension of the work array workl
        lworkl_ = ncv_ * (3 * ncv_ + 5);
    } else {
        std::cout<<"PARPACK Solver only defined for double & cdouble!\n";
    }
    // default tol: machine precision
    tol_ = 0.0;
    sigma_ = 0.0;

    // work array init to 0
    workd_pt = std::vector<T>(3 * nlocmax_, 0.0);
    // work array of length lworkl init to 0
    workl_pt = std::vector<T>(lworkl_, 0.0);
    V_pt = std::vector<T>(ncv_ * nlocmax_, 0.0);
    d_pt = std::vector<T>(nev_ + 1, 0.0);
    z_pt = std::vector<T>((nlocmax_ + 1) * (nev_ + 1), 0.0);
    // residual vector. info = 0 random initial vector is used or resid will be used.
    resid_pt = std::vector<T>(nlocmax_, 0.0);
    select_pt = std::vector<a_int>(ncv_, 1);
    
    if constexpr (std::is_same<T, cdouble>::value) {
        rwork_pt = std::vector<T>(ncv_, 0.0);
        workev_pt = std::vector<T>(2*ncv_, 0.0);
    }
    
    // method for selecting the implicit shifts
    iparam_[0] = 1;
    // On Input: maximum number of Arnoldi update iterations allowed;
    // On output: actual number of interations taken
    iparam_[2] = PARPACK_MAXITERATION;
    // only work for 1
    iparam_[3] = 1;
    // number of ev found by arpack. number of ev satisfy convengence criterion
    iparam_[4] = 0;
    // iparam_[4] = ncv_;
    // Mode 1: A*x=lambda*x
    iparam_[6] = 1;
}

template <class T>
void PARPACKSolver<T>::setStartVec(T* vec) {
    info_ = 1;
    #pragma omp parallel for
    for (idx_t i = 0; i < nloc_; ++i) {
        resid_pt[i] = vec[i];
    }
}

template <class T>
void PARPACKSolver<T>::setNev(a_int nev) {
    nev_ = nev;
}
template <class T>
void PARPACKSolver<T>::setNcv(a_int ncv) {
    ncv_ = ncv;
}
template <class T>
void PARPACKSolver<T>::setMaxIter(int iterNum) {
    iparam_[2] = iterNum;
}

template <class T>
void PARPACKSolver<T>::setRvec(a_int rvec) {
    rvec_ = rvec;
}

template <class T>
void PARPACKSolver<T>::check( ) {
    assert_msg(M_, "Matrix is not set for PARPACK Solver!");
    assert_msg(nev_ > 0, "nev should be a positive integer for PARPACK Solver!");
}

template <class T>
void PARPACKSolver<T>::diag( ) {
    check();
    Timer timer;
    int workerID = M_->getWorkerID();
    bool ismaster = (workerID == MPI_MASTER);
    if (ismaster) std::cout<<"Begin PARPACK Iteration and timer started...\n";
    timer.tik();
    run();
    timer.tok();
    if (ismaster) {
        std::cout<<"INFO:"<<info_<<". Total iteration:"<<iparam_[2]<<". Total time:"<<timer.elapse()<<" milliseconds.\n";
    }

    if (ismaster) std::cout<<"Begin post processing...\n";
    timer.tik();
    postRun();
    timer.tok();
    if (ismaster) {
        std::cout<<"INFO:"<<info_<<". Total ev found:"<<iparam_[4]<<". post processing time:"<<timer.elapse()<<" milliseconds.\n";
    }
    if (ismaster) {
        std::cout<<"Eigenvalues: ";
        for (int i = 0; i < nev_; ++i) std::cout<<d_pt[i]<<", ";
        std::cout<<"\n";
    }
}

template <class T>
void PARPACKSolver<T>::diag(double spin, SSOp<T>* SS){
    check();
    Timer timer;
    int workerID = M_->getWorkerID();
    bool ismaster = (workerID == MPI_MASTER);
    if (ismaster) std::cout<<"Begin PARPACK Iteration and timer started...\n";
    timer.tik();
    run(spin, SS);
    timer.tok();
    if (ismaster) {
        std::cout<<"INFO:"<<info_<<". Total iteration:"<<iparam_[2]<<". Total time:"<<timer.elapse()<<" milliseconds\n";
    }
    if (ismaster) std::cout<<"Begin post processing..."<<std::endl;
    timer.tik();
    postRun();
    timer.tok();
    if (ismaster) {
        std::cout<<"INFO:"<<info_<<". Total ev found:"<<iparam_[4]<<". post processing time:"<<timer.elapse()<<" milliseconds\n";
    }
    if (ismaster) {
        std::cout<<"Eigenvalues: ";
        for (int i = 0; i < nev_; ++i) std::cout<<d_pt[i]<<", ";
        std::cout<<"\n";
    }
}

template <class T>
void PARPACKSolver<T>::run( ) {
    while (ido_ != 99) {
        if constexpr (std::is_same<T, double>::value) {
            arpack::saupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_algebraic, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, info_);
        } else if constexpr (std::is_same<T, cdouble>::value) { 
            arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_realpart, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, rwork_pt.data(), info_);
        }
        M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
    }

    // check number of ev found by arpack
    if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
        std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<"\n";
        exit(1);
    }
}

template <class T>
void PARPACKSolver<T>::run(double spin, SSOp<T>* SS){
    while (ido_ != 99) {
        if constexpr (std::is_same<T, double>::value) {
            arpack::saupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_algebraic, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, info_);
        } else if constexpr (std::is_same<T, cdouble>::value) { 
            arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_realpart, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, rwork_pt.data(), info_);
        }
        SS->project(spin, &(workd_pt[ipntr_[0] - 1]));
        M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
    }
    // check number of ev found by arpack
    if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
        std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ <<",iterations taken:"<<iparam_[2]<< std::endl;
        exit(1);
    }
}

// project out solved eigen states
template <class T>
void PARPACKSolver<T>::run(T* states, int statesNum, double penalty){
    while (ido_ != 99) {
        if constexpr (std::is_same<T, double>::value) {
            arpack::saupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_algebraic, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, info_);
        } else if constexpr (std::is_same<T, cdouble>::value) { 
            arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_, arpack::which::smallest_realpart, nev_, tol_, resid_pt.data(), ncv_, \
                V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, rwork_pt.data(), info_);
        }
        M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
        for (int i = 0; i < statesNum; ++i){
            cdouble overlap;
            vConjDotv(states+i*nloc_, &(workd_pt[ipntr_[0]-1], overlap, nloc_));
            saxpy(&(workd_pt[ipntr_[1] - 1]), overlap*penalty, states+i*nloc_, nloc_);
        }
    }

  // check number of ev found by arpack
  if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
    std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
    exit(1);
  }
}

template <class T>
void PARPACKSolver<T>::postRun( ) {
    if constexpr (std::is_same<T, double>::value) {
        arpack::seupd(MCW, rvec_, arpack::howmny::ritz_vectors, select_pt.data(), d_pt.data(), z_pt.data(), ldz_, sigma_, arpack::bmat::identity, nloc_, \
            arpack::which::smallest_algebraic, nev_, tol_, resid_pt.data(), ncv_, V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, info_);
    } else if constexpr (std::is_same<T, cdouble>::value) {
        arpack::neupd(MCW, rvec_, arpack::howmny::ritz_vectors, select_pt.data(), d_pt.data(), z_pt.data(), ldz_, sigma_, workev_pt.data(), arpack::bmat::identity, nloc_, \
            arpack::which::smallest_realpart, nev_, tol_, resid_pt.data(), ncv_, V_pt.data(), ldv_, iparam_.data(), ipntr_.data(), workd_pt.data(), workl_pt.data(), lworkl_, rwork_pt.data(), info_);
    }
    if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
    std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
    exit(1);
  }
}

// /*
//     *********************************************
//     * PARPACK Solver for real symmetric problem *
//     *********************************************
// */
// // T is float or double
// template <class T>
// class PARPACKRealSolver{
// public:
//     MPI_Fint MCW;
// //    SparseMatrix<T> *M_; // pointer to the SparseMatrix object
//     BaseMatrix<T> *M_;
//     a_uint ntot_, nlocmax_, nloc_, ldz_, ldv_;
//     int rank_;
//     a_int ido_, info_, nev_, ncv_, lworkl_, rvec_;
//     T tol_, sigma_;
//     std::array<a_int, 14> ipntr_; // output. ipntr(0) pointer to current operand vector in workd. ipntr(1) pointer to current result vector in workd
//     std::array<a_int, 11> iparam_;
    
//     T *workd_pt, *workl_pt, *V_pt, *d_pt, *z_pt, *resid_pt;
//     a_int *select_pt;
    

//     PARPACKRealSolver(){};
//     PARPACKRealSolver(BaseMatrix<T> *M, a_int nev);
//     ~PARPACKRealSolver();
    
//     T* getEigval() const {return d_pt;}
//     T* getEigvec(int ithVec=0) const {return V_pt+nloc_*ithVec;}
//     void setNev(a_int nev);
//     void setNcv(a_int ncv);
//     void setMaxIter(int iterNum);
//     void setRvec(a_int rvec);
        
//     void run();
//     void postRun();  
//     void diag();   
// };

// /*
//     ****************************************
//     * PARPACK Solver for Hermitian problem *
//     ****************************************
// */
// // T is float or double
// template <class T>
// class PARPACKComplexSolver{
//     MPI_Fint MCW;
// //        SparseMatrix<std::complex<T>> *M_; // pointer to the SparseMatrix object
//     BaseMatrix<std::complex<T>> *M_;
//     a_uint ntot_, nlocmax_, nloc_, ldz_, ldv_;
//     int rank_;
//     a_int ido_, info_, nev_, ncv_, lworkl_, rvec_;
//     T tol_;
//     std::complex<T> sigma_;
//     std::array<a_int, 14> ipntr_; // output. ipntr(0) pointer to current operand vector in workd. ipntr(1) pointer to current result vector in workd
//     std::array<a_int, 11> iparam_;
// `
//     std::complex<T> *workd_pt, *workl_pt, *V_pt, *d_pt, *z_pt, *resid_pt;
//     std::complex<T> *rwork_pt, *workev_pt;
//     a_int *select_pt;

// public:
//     PARPACKComplexSolver(){};
//     PARPACKComplexSolver(BaseMatrix<std::complex<T>> *M, a_int nev);
//     ~PARPACKComplexSolver();

//     std::complex<T>* getEigval() const {return d_pt;}
//     std::complex<T>* getEigvec(int ithVec=0) const {return V_pt+nloc_*ithVec;}
//     void setStartVec(std::complex<T>* vec);
//     void setNev(a_int nev);
//     void setNcv(a_int ncv);
//     void setMaxIter(int iterNum);
//     void setRvec(a_int rvec);

//     void run();
//     void run(double spin,SSOp<std::complex<T>>* SS);
//     void run(std::complex<T>* states, int statesNum, double penalty=1000.0);
//     void postRun();
//     void diag();
//     void diag(double spin,SSOp<std::complex<T>>* SS);
// };

// /*
//     ******************
//     * Real Symmetric *
//     ******************
// */
// template <class T>
// PARPACKRealSolver<T>::PARPACKRealSolver(BaseMatrix<T> *H, a_int nev){
//     MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
//     info_ = 0; //Input:use random initial vector. Output:0-normal exit,1-maximum iteration taken
//     ido_ = 0; // reverse communication flag 
    
//     M_ = H; // pointer to the hamitonian object
//     ntot_ = M_->getntot();
//     nlocmax_ = M_->getnlocmax(); // maximum local vector size
//     nloc_ = M_->getnloc(); // actual local vector size
//     nev_ = nev; // number of eigenvalues to be computed
//     ncv_ = (2 * nev_ + 1) > PARPACK_MINNCV?(2 * nev_ + 1):PARPACK_MINNCV; // number of Lanczos vectors generated at each iteration. number of colunm of V
//     if (ncv_ > M_->getDim()) ncv_ = M_->getDim()-1;
//     ldz_ = nlocmax_ + 1;
//     lworkl_ = ncv_ * (ncv_ + 8); // dimension of the work array workl
//     ldv_ = nlocmax_; // leading dimension of V
//     rvec_ = 1;
    
//     tol_ = (T) 0.0; // machine precision
//     sigma_ = (T) 0.0;

//     workd_pt = new(std::nothrow) T[3 * nlocmax_]; assert(workd_pt!=NULL);// work array initial to 0
//     workl_pt = new(std::nothrow) T[lworkl_]; assert(workl_pt!=NULL);// work array of length lworkl initial to 0
//     V_pt = new(std::nothrow) T[ncv_ * nlocmax_]; assert(V_pt!=NULL);//output:store eigenvectors
//     d_pt = new(std::nothrow) T[nev_ + 1]; assert(d_pt!=NULL);//output:store eigenvalues
//     z_pt = new(std::nothrow) T[(nlocmax_ + 1) * (nev_ + 1)]; assert(z_pt!=NULL);
//     resid_pt = new(std::nothrow) T[nlocmax_]; assert(resid_pt!=NULL);// residual vector. info=0 random initial vector is used or resid will be used.
//     select_pt = new(std::nothrow) a_int[ncv_]; assert(select_pt!=NULL);
    
// //     vecInit<T, a_uint>(workd_pt, 3 * nlocmax_, 0.0);
// //     vecInit<T, a_int>(workl_pt, lworkl_, 0.0);
// //     vecInit<a_int, a_int>(ipntr_.data(), 14, 0);
// //     vecInit<a_int, a_int>(iparam_.data(), 11, 0);
//     vecInit<a_int, a_int>(select_pt, ncv_, 1);

//     iparam_[0] = 1; // method for selecting the implicit shifts
//     iparam_[2] = PARPACK_MAXITERATION; // On Input: maximum number of Arnoldi update iterations allowed; On output: actual number of interations taken
//     iparam_[3] = 1; // only work for 1
//     iparam_[4] = 0;  // number of ev found by arpack. number of ev satisfy convengence criterion
//     iparam_[6] = 1; // Mode 1: A*x=lambda*x
    
// }

// template <class T>
// PARPACKRealSolver<T>::~PARPACKRealSolver(){
//     delete [] workd_pt;
//     delete [] workl_pt;
//     delete [] V_pt;
//     delete [] d_pt;
//     delete [] z_pt;
//     delete [] resid_pt;
//     delete [] select_pt;
// }

// template <class T>
// void PARPACKRealSolver<T>::run(){
//     while (ido_ != 99) {
//         arpack::saupd(MCW, ido_, arpack::bmat::identity, nloc_,
//                     arpack::which::smallest_algebraic, nev_, tol_, resid_pt, ncv_,
//                     V_pt, ldv_, iparam_.data(), ipntr_.data(), workd_pt,
//                     workl_pt, lworkl_, info_);
//         M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
//     }
//     // check number of ev found by arpack
//     if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//         std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ <<",iterations taken:"<<iparam_[2]<< std::endl;
//         exit(1);
//     }
// }

// template <class T>
// void PARPACKRealSolver<T>::postRun(){
//     arpack::seupd(MCW, rvec_, arpack::howmny::ritz_vectors, select_pt,
//                 d_pt, z_pt, ldz_, sigma_, arpack::bmat::identity, nloc_,
//                 arpack::which::smallest_algebraic, nev_, tol_, resid_pt, ncv_,
//                 V_pt, ldv_, iparam_.data(), ipntr_.data(), workd_pt, workl_pt, lworkl_, info_);
//     if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//     std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
//     exit(1);
//   }
// }

// template <class T>
// void PARPACKRealSolver<T>::diag(){
//     Timer timer;
//     int workerID = M_->getWorkerID();
//     if (workerID==MPI_MASTER) std::cout<<"Begin PARPACK Iteration and timer started..."<<std::endl;
//     timer.tik();
//     run();
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total iteration:"<<iparam_[2]<<". Total time:"<<timer.elapse()<<" milliseconds"<<std::endl;

//     if (workerID==MPI_MASTER) std::cout<<"Begin post processing..."<<std::endl;
//     timer.tik();
//     postRun();
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total ev found:"<<iparam_[4]<<". post processing time:"<<timer.elapse()<<" milliseconds"<<std::endl;
//     if (workerID==MPI_MASTER) {
//         std::cout<<"Eigenvalues: ";
//         for (int i = 0; i < nev_; ++i) std::cout<<d_pt[i]<<", ";
//         std::cout<<std::endl;
//     }
// }

// template <class T>
// void PARPACKRealSolver<T>::setNev(a_int nev){
//     nev_ = nev;
// }
    
// template <class T>
// void PARPACKRealSolver<T>::setNcv(a_int ncv){
//     ncv_ = ncv;
// }

// template <class T>
// void PARPACKRealSolver<T>::setMaxIter(int iterNum){
//     iparam_[2] = iterNum;
// }

// template <class T>
// void PARPACKRealSolver<T>::setRvec(a_int rvec){
//     rvec_ = rvec;
// }

// /*
//     *************
//     * Hermitian *
//     *************
// */
// template <class T>
// PARPACKComplexSolver<T>::PARPACKComplexSolver(BaseMatrix<std::complex<T>> *H, a_int nev){
//     if (H->getDim() <= 0) {
//         std::cout<<"Empty Matrix for PARPACK Solver!\n";
//         exit(1);
//     }
//     MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
//     info_ = 0; //Input:use random initial vector. Output:0-normal exit,1-maximum iteration taken
//     ido_ = 0; // reverse communication flag 
    
//     M_ = H; // pointer to the hamitonian object
//     ntot_ = M_->getntot();
//     nlocmax_ = M_->getnlocmax(); // maximum local vector size
//     nloc_ = M_->getnloc(); // actual local vector size
//     nev_ = nev; // number of eigenvalues to be computed
//     ncv_ = (2 * nev_ + 1) > PARPACK_MINNCV?(2 * nev_ + 1):PARPACK_MINNCV; // number of Lanczos vectors generated at each iteration. number of colunm of V
//     assert_msg(nev_ > 0, "nev_ <= 0");
//     assert_msg(ncv_ > 0, "ncv_ <= 0");
//     if (ncv_ > M_->getDim()) ncv_ = M_->getDim()-1;
//     ldz_ = nlocmax_ + 1;
//     lworkl_ = ncv_ * (3 * ncv_ + 5); // dimension of the work array workl
//     ldv_ = nlocmax_; // leading dimension of V
//     rvec_ = 1;
    
//     tol_ = (T) 0.0; // machine precision
//     sigma_ = {0.0, 0.0};

//     workd_pt = new(std::nothrow) std::complex<T>[3 * nlocmax_]; assert(workd_pt!=NULL);// work array initial to 0
//     workl_pt = new(std::nothrow) std::complex<T>[lworkl_]; assert(workl_pt!=NULL);// work array of length lworkl initial to 0
//     V_pt = new(std::nothrow) std::complex<T>[ncv_ * nlocmax_]; assert(V_pt!=NULL);
//     d_pt = new(std::nothrow) std::complex<T>[nev_ + 1]; assert(d_pt!=NULL);
//     z_pt = new(std::nothrow) std::complex<T>[(nlocmax_ + 1) * (nev_ + 1)];assert(z_pt!=NULL);
//     resid_pt = new(std::nothrow) std::complex<T>[nlocmax_]; assert(resid_pt!=NULL);// residual vector. info=0 random initial vector is used or resid will be used.
//     select_pt = new(std::nothrow) a_int[ncv_]; assert(select_pt!=NULL);
    
//     rwork_pt = new(std::nothrow) std::complex<T>[ncv_]; assert(rwork_pt!=NULL);
//     workev_pt = new(std::nothrow) std::complex<T>[2*ncv_]; assert(workev_pt!=NULL);
    
// //     vecInit<T, a_uint>(workd_pt, 3 * nlocmax_, 0.0);
// //     vecInit<T, a_int>(workl_pt, lworkl_, 0.0);
// //     vecInit<a_int, a_int>(ipntr_.data(), 14, 0);
// //     vecInit<a_int, a_int>(iparam_.data(), 11, 0);
//     vecInit<a_int, a_int>(select_pt, ncv_, 1);

//     iparam_[0] = 1; // method for selecting the implicit shifts
//     iparam_[2] = PARPACK_MAXITERATION; // On Input: maximum number of Arnoldi update iterations allowed; On output: actual number of interations taken
//     iparam_[3] = 1; // only work for 1
//     iparam_[4] = 0;  // number of ev found by arpack. number of ev satisfy convengence criterion
// //     iparam_[4] = ncv_;
//     iparam_[6] = 1; // Mode 1: A*x=lambda*x
    
// }

// template <class T>
// void PARPACKComplexSolver<T>::setStartVec(std::complex<T>* vec){
//     info_ = 1;
//     #pragma omp parallel for
//     for (idx_t i = 0; i < nloc_; ++i){
//         resid_pt[i] = vec[i];
//     }
// }

// template <class T>
// PARPACKComplexSolver<T>::~PARPACKComplexSolver(){
//     delete [] workd_pt;
//     delete [] workl_pt;
//     delete [] V_pt;
//     delete [] d_pt;
//     delete [] z_pt;
//     delete [] resid_pt;
//     delete [] select_pt;
//     delete [] rwork_pt;
//     delete [] workev_pt;
// }

// template <class T>
// void PARPACKComplexSolver<T>::run(){
//     while (ido_ != 99) {
//     arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_,
//                   arpack::which::smallest_realpart, nev_, tol_, resid_pt, ncv_,
//                   V_pt, ldv_, iparam_.data(), ipntr_.data(), workd_pt,
//                   workl_pt, lworkl_, rwork_pt, info_);
//     M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
//   }

//   // check number of ev found by arpack
//   if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//     std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
//     exit(1);
//   }
// }
// template <class T>
// void PARPACKComplexSolver<T>::run(double spin,SSOp<std::complex<T>>* SS){
//     while (ido_ != 99) {
//         arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_,
//                   arpack::which::smallest_realpart, nev_, tol_, resid_pt, ncv_,
//                   V_pt, ldv_, iparam_.data(), ipntr_.data(), workd_pt,
//                   workl_pt, lworkl_, rwork_pt, info_);
//         SS->project(spin, &(workd_pt[ipntr_[0] - 1]));
//         M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
//     }
//     // check number of ev found by arpack
//     if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//         std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ <<",iterations taken:"<<iparam_[2]<< std::endl;
//         exit(1);
//     }
// }

// // project out solved eigen states
// template <class T>
// void PARPACKComplexSolver<T>::run(std::complex<T>* states, int statesNum, double penalty){
//     while (ido_ != 99) {
//     arpack::naupd(MCW, ido_, arpack::bmat::identity, nloc_,
//                   arpack::which::smallest_realpart, nev_, tol_, resid_pt, ncv_,
//                   V_pt, ldv_, iparam_.data(), ipntr_.data(), workd_pt,
//                   workl_pt, lworkl_, rwork_pt, info_);
//     M_->MxV(&(workd_pt[ipntr_[0] - 1]), &(workd_pt[ipntr_[1] - 1]));
//     for (int i = 0; i < statesNum; ++i){
//         cdouble overlap;
//         vConjDotv(states+i*nloc_, &(workd_pt[ipntr_[0]-1], overlap, nloc_));
//         saxpy(&(workd_pt[ipntr_[1] - 1]), overlap*penalty, states+i*nloc_, nloc_);
//     }
//   }

//   // check number of ev found by arpack
//   if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//     std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
//     exit(1);
//   }
// }

// template <class T>
// void PARPACKComplexSolver<T>::postRun(){
//     arpack::neupd(MCW, rvec_, arpack::howmny::ritz_vectors, select_pt,
//                 d_pt, z_pt, ldz_, sigma_, workev_pt,
//                 arpack::bmat::identity, nloc_, arpack::which::smallest_realpart,
//                 nev_, tol_, resid_pt, ncv_, V_pt, ldv_, iparam_.data(),
//                 ipntr_.data(), workd_pt, workl_pt, lworkl_, rwork_pt,
//                 info_);
//     if (iparam_[4] < nev_ /*arpack may succeed to compute more EV than expected*/ || info_ != 0) {
//     std::cout << "ERROR: iparam[4] " << iparam_[4] << ", nev " << nev_ << ", info " << info_ << ",iterations taken:"<<iparam_[2]<<std::endl;
//     exit(1);
//   }
// }

// template <class T>
// void PARPACKComplexSolver<T>::setNev(a_int nev){
//     nev_ = nev;
// }
// template <class T>
// void PARPACKComplexSolver<T>::setNcv(a_int ncv){
//     ncv_ = ncv;
// }
// template <class T>
// void PARPACKComplexSolver<T>::setMaxIter(int iterNum){
//     iparam_[2] = iterNum;
// }

// template <class T>
// void PARPACKComplexSolver<T>::setRvec(a_int rvec){
//     rvec_ = rvec;
// }

// template <class T>
// void PARPACKComplexSolver<T>::diag(){
//     Timer timer;
//     int workerID = M_->getWorkerID();
//     if (workerID==MPI_MASTER) std::cout<<"Begin PARPACK Iteration and timer started..."<<std::endl;
//     timer.tik();
//     run();
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total iteration:"<<iparam_[2]<<". Total time:"<<timer.elapse()<<" milliseconds"<<std::endl;

//     if (workerID==MPI_MASTER) std::cout<<"Begin post processing..."<<std::endl;
//     timer.tik();
//     postRun();
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total ev found:"<<iparam_[4]<<". post processing time:"<<timer.elapse()<<" milliseconds"<<std::endl;
//     if (workerID==MPI_MASTER) {
//         std::cout<<"Eigenvalues: ";
//         for (int i = 0; i < nev_; ++i) std::cout<<d_pt[i]<<", ";
//         std::cout<<std::endl;
//     }
// }

// template <class T>
// void PARPACKComplexSolver<T>::diag(double spin, SSOp<std::complex<T>>* SS){
//     Timer timer;
//     int workerID = M_->getWorkerID();
//     if (workerID==MPI_MASTER) std::cout<<"Begin PARPACK Iteration and timer started..."<<std::endl;
//     timer.tik();
//     run(spin, SS);
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total iteration:"<<iparam_[2]<<". Total time:"<<timer.elapse()<<" milliseconds"<<std::endl;

//     if (workerID==MPI_MASTER) std::cout<<"Begin post processing..."<<std::endl;
//     timer.tik();
//     postRun();
//     timer.tok();
//     if (workerID==MPI_MASTER) std::cout<<"INFO:"<<info_<<". Total ev found:"<<iparam_[4]<<". post processing time:"<<timer.elapse()<<" milliseconds"<<std::endl;
//     if (workerID==MPI_MASTER) {
//         std::cout<<"Eigenvalues: ";
//         for (int i = 0; i < nev_; ++i) std::cout<<d_pt[i]<<", ";
//         std::cout<<std::endl;
//     }
// }

#endif // PARPACKSolver_hpp
