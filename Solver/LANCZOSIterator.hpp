//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
//

#ifndef LANCZOSIterator_hpp
#define LANCZOSIterator_hpp

#include <vector>
#include <cmath>

#include "Global/globalType.hpp"
#include "Operator/SparseMatrix.hpp"
#include "Algebra/algebra.hpp"

// // #define USE_MKL

// #ifdef USE_MKL

// template <class T>
// class LANCZOSIterator{
// public:
//     LANCZOSIterator( ) { }
//     LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt = LANCZOS_DEFAULT_OPTION);
//     ~LANCZOSIterator( ) { }
    
//     void init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt);
//     void run(T* vec);
//     // reorthogonalize
//     void runReOth(T* vec);

// protected:
//     // local vectorsize
//     BaseMatrix<T> *M_{nullptr};
//     LANCZOS_OPTION option{LANCZOS_OPTION::ALPHA_BETA};
//     int krylovDim;
//     // stores alpha, beta
//     std::vector<double> alpha, beta;
//     double vecNorm;
//     // q_j-1, q_j, r
//     std::vector<T> q_pre, q, r;
//     // AQ = QH + residual
//     std::vector<T> Q;
// };

// template <class T>
// LANCZOSIterator<T>::LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
//     init(M, krydim, opt);
// }

// template <class T>
// void LANCZOSIterator<T>::init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
//     assert_msg(M->getDim()== M->getColDim(), "M must be a square matrix in LANCZOSIterator!");
//     M_ = M;
//     option = opt;
//     krylovDim = krydim;
//     q_pre = std::vector<T>(M->getnlocmax(), 0.0);
//     q = std::vector<T>(M->getnlocmax(), 0.0);
//     r = std::vector<T>(M->getnlocmax(), 0.0);
//     if (option == ALPHA_BETA_Q) {
//         Q = std::vector<T>(krydim * M->getnloc());
//     }
// }

// template <class T>
// void LANCZOSIterator<T>::run(T* vec){
//     alpha.clear();
//     beta.clear();
//     T result;
//     T zero{0.0}, one{1.0};
//     auto nloc = M_->getnloc(); 
//     // q = vec/|vec|
//     MKL::copy(nloc, vec, q.data());
//     vecNorm = MKL::mpiNorm(nloc, vec);
//     if (vecNorm < INFINITESIMAL){
//         if (isMaster) std::cout<<"norm(vecIn) = "<<vecNorm<<". Check if the input vec is 0.\n";
//         return;
//     }
//     MKL::scale(nloc, 1.0/vecNorm, q.data());
//     // Q1 = [q]
//     if (option == ALPHA_BETA_Q) MKL::copy(nloc, q.data(), Q.data());
//     // r = Aq
//     M_->MxV(q.data(), r.data());
//     // a1 = q^r
//     alpha.push_back(std::real(MKL::mpiDot(nloc, q.data(), r.data())));
//     // r = r - a1*q
//     MKL::axpy(nloc, -(T)alpha[0], q.data(), r.data());
//     // beta[1] = |r|
//     beta.push_back(MKL::mpiNorm(nloc, r.data()));
//     bool isMaster = (M_->getWorkerID() == MPI_MASTER);
//     if ( isMaster && beta[0] == 0.0){
//         std::cout<<"LANZOSIteration stops for beta[0] = 0.\n";
//         return;
//     }
    
//     for (int i = 1; i < krylovDim; ++i){
//         // q_pre = q
//         MKL::copy(nloc, q.data(), q_pre.data());
//         // q = r/beta[i-1]
//         MKL::scale(nloc, one/beta[i-1], r.data());
//         MKL::copy(nloc, r.data(), q.data());
//         //Qi = [Qi-1,q]
//         if (option == ALPHA_BETA_Q) MKL::copy(nloc, q.data(), Q.data() + i * nloc);
//         // r = Aq
//         M_->MxV(q.data(), r.data());
//         //r = r - beta[i-1]*q_pre
//         MKL::axpy(nloc, -(T)beta[i-1], q_pre.data(), r.data());
//         //alpha[i] = q^r
//         alpha.push_back(std::real(MKL::mpiDot(nloc, q.data(), r.data())));
//         // r = r - alpha[i]*q
//         MKL::axpy(nloc, -(T)alpha[i], q.data(), r.data());
//         // beta[i] = |r|
//         beta.push_back(MKL::mpiNorm(nloc, r.data()));
//         if (isMaster && beta[i] == zero){
//             std::cout<<"LANZOSIteration stops for beta["<<i<<"] = 0.\n";
//             return;
//         }
//     }
// }

// #else

template <class T>
class LANCZOSIterator{
public:
    LANCZOSIterator( ) { }
    LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt = LANCZOS_OPTION::ALPHA_BETA);
    ~LANCZOSIterator( );
    
    void init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt);
    void run(T* vec);
    // reorthogonalize
    void runReOth(T* vec);

protected:
    // local vectorsize
    BaseMatrix<T> *M_{nullptr};
    LANCZOS_OPTION option{LANCZOS_OPTION::ALPHA_BETA};
    int krylovDim;
    // stores alpha, beta
    std::vector<double> alpha, beta;
    double vecNorm;
    // q_j-1, q_j, r
    cdouble *q_pre{nullptr}, *q{nullptr}, *r{nullptr};
    // AQ = QH + residual
    cdouble *Q{nullptr};
};

template <class T>
LANCZOSIterator<T>::LANCZOSIterator(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
    init(M, krydim, opt);
}

template <class T>
void LANCZOSIterator<T>::init(BaseMatrix<T> *M, int krydim, LANCZOS_OPTION opt){
    assert_msg(M->getDim()== M->getColDim(), "M must be a square matrix in LANCZOSIterator!");
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
void LANCZOSIterator<T>::run(T* vec){
    alpha.clear();
    beta.clear();
    auto nloc =  M_->getnloc();
    bool isMaster = (M_->getWorkerID() == MPI_MASTER);
    // q = vec/|vec|
    T zero{0.0}, one{1.0};
    copy(nloc, vec, q);
    vecNorm = mpiNorm(vec, nloc);
    if (vecNorm < INFINITESIMAL){
        if (isMaster) std::cout<<"norm(vecIn) = "<<vecNorm<<". Check if the input vec is 0.\n";
        return;
    }
    scale(q, 1.0/(T)vecNorm,nloc);
    // Q1 = [q]
    if (option == ALPHA_BETA_Q) copy(nloc, q, &Q[0]);
    // r = Aq
    M_->MxV(q, r);
    // a1 = q^r
    alpha.push_back(std::real(mpiDot(q, r, nloc)));
    // r = r - a1*q
    axpy(r, -(cdouble)alpha[0], q, nloc);
    // beta[1] = |r|
    beta.push_back(mpiNorm(r, nloc));
    if (isMaster && beta[0]==0.0){
        std::cout<<"LANZOSIteration stops for beta[0] = 0.\n";
        return;
    }
    
    for (int i = 1; i < krylovDim; ++i){
        // q_pre = q
        copy(nloc, q, q_pre);
        // q = r/beta[i-1]
        copy(q, one/beta[i-1], r,nloc);
        //Qi = [Qi-1,q]
        if (option == ALPHA_BETA_Q) copy(nloc, q, &Q[i * nloc]);
        // r = Aq
        M_->MxV(q, r);
        //r = r - beta[i-1]*q_pre
        axpy(r, -(cdouble)beta[i-1], q_pre, nloc);
        //alpha[i] = q^r
        alpha.push_back(std::real(mpiDot(q, r, nloc)));
        // r = r - alpha[i]*q
        axpy(r, -(cdouble)alpha[i], q, nloc);
        // beta[i] = |r|
        beta.push_back(mpiNorm(r, nloc));
        if (isMaster && beta[i] == zero){
            std::cout<<"LANZOSIteration stops for beta["<<i<<"] = 0."<<std::endl;
            return;
        }
    }
}

// #endif // USE_MKL

#endif // LANCZOSIterator_hpp
