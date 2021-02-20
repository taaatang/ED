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

    // BLAS 1
    // y = a * x + y
    void axpy(MKL_INT n, const double a, const double *x, double *y) {
        cblas_daxpy(n, a, x, 1, y, 1);
    }
    void axpy(const MKL_INT n, const cdouble a, const cdouble *x, cdouble *y) {
        cblas_zaxpy(n, a, x, 1, y,1);
    }
    // y = x
    void copy(const MKL_INT n, const double *x, double *y) {
        cblas_dcopy(n, x, 1, y, 1);
    }
    void copy(const MKL_INT n, const cdouble *x, cdouble *y) {
        cblas_zcopy(n, x, 1, y, 1);
    }
    // conj(x).dot(y)
    double dot(const MKL_INT n, const double *x, const double *y) {
        cblas_ddot(n, x, 1, y, 1);
    }
    cdouble dot(const MKL_INT n, const cdouble *x, const cdouble *y) {
        cdouble res = 0.0;
        cblas_zdotc_sub(n, x, 1, y, q, &res);
        return res;
    }
    // Euclidean norm
    double norm(const MKL_INT n, const double *x){
        cblas_dnrm2(n, x, 1);
    }
    double norm(const MKL_INT n, const cdouble *x){
        cblas_dznrm2(n, x, 1);
    }
    // x = a * x
    void scale(const MKL_INT n, const double a, double *x) {
        cblas_dscal(n, a, x, 1);
    }
    void scale(const MKL_INT n, const cdouble a, cdouble *x) {
        cblas_zscal(n, a, x, 1);
    }
    void scale(const MKL_INT n, const double a, cdouble *x) {
        cblas_zdscal(n, a, x, 1);
    }
}


/*
    Eigen solver for a self adjoint matrix
*/
#endif
