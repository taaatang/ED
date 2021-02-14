//
// algebra.hpp
// ED
//
// Created by tatang on 1/15/2020
//
//

/*
    C++ interface to mkl library
*/
#ifndef algebra_hpp
#define algebra_hpp
// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>
#include "Global/globalType.hpp"
#define _MKL_
#define MKL_INT idx_t
#define MKL_Complex16 cdouble
#include "mkl.h"
#include "omp.h"
#include <vector>

namespace MKL{
    /*
        LARPACK eigen solver for a real symmetric tridiagonal matrix
    */
    void diagTri(std::vector<double>* a, std::vector<double>* b, std::vector<double>* U);
    /*
        inspector-executer mkl sparse.
    */
    void create(sparse_matrix_t& A, idx_t rows, idx_t cols, std::vector<MKL_INT>& rowInitList, std::vector<MKL_INT>& colList, std::vector<cdouble>& valList, MKL_INT mvNum=100);

    void MxV(sparse_matrix_t& A, MKL_Complex16* vin, MKL_Complex16* vout, MKL_Complex16 alpha=1.0, MKL_Complex16 beta=1.0);

    void create(sparse_matrix_t& A, idx_t rows, idx_t cols, std::vector<MKL_INT>& rowInitList, std::vector<MKL_INT>& colList, std::vector<double>& valList, MKL_INT mvNum=100);
    void MxV(sparse_matrix_t& A, double* vin, double* vout, double alpha=1.0, double beta=1.0);

    void destroy(sparse_matrix_t& A);
}


/*
    Eigen solver for a self adjoint matrix
*/
#endif
