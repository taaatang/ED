//
// algebra.hpp
// ED
//
// Created by tatang on 1/15/2020
//
//

#ifndef algebra_hpp
#define algebra_hpp

// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include "omp.h"

#include "Global/globalType.hpp"
#define MKL_INT idx_t
#define MKL_Complex16 cdouble
#include "mkl.h"

#include "Utils/mpiwrap.hpp"
#include "Utils/io.hpp"

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
    inline void axpy(MKL_INT n, const double a, const double *x, double *y) {
        cblas_daxpy(n, a, x, 1, y, 1);
    }
    inline void axpy(const MKL_INT n, const cdouble a, const cdouble *x, cdouble *y) {
        cblas_zaxpy(n, &a, x, 1, y,1);
    }
    // y = x
    inline void copy(const MKL_INT n, const double *x, double *y) {
        cblas_dcopy(n, x, 1, y, 1);
    }
    inline void copy(const MKL_INT n, const cdouble *x, cdouble *y) {
        cblas_zcopy(n, x, 1, y, 1);
    }
    // conj(x).dot(y)
    inline double dot(const MKL_INT n, const double *x, const double *y) {
        return cblas_ddot(n, x, 1, y, 1);
    }
    inline cdouble dot(const MKL_INT n, const cdouble *x, const cdouble *y) {
        cdouble res = 0.0;
        cblas_zdotc_sub(n, x, 1, y, 1, &res);
        return res;
    }
    template <typename T>
    T mpiDot(const MKL_INT n, const T *x, const T *y) {
        T resLoc = dot(n, x, y);
        T res;
        MPI_Allreduce(&resLoc, &res, 1);
        return res;
    }
    // Euclidean norm
    inline double norm(const MKL_INT n, const double *x) {
        return cblas_dnrm2(n, x, 1);
    }
    inline double norm(const MKL_INT n, const cdouble *x) {
        return cblas_dznrm2(n, x, 1);
    }
    template <typename T>
    double mpiNorm(const MKL_INT n, const T *x) {
        return std::real(std::sqrt(mpiDot(n, x, x)));
    }
    // x = a * x
    inline void scale(const MKL_INT n, const double a, double *x) {
        cblas_dscal(n, a, x, 1);
    }
    inline void scale(const MKL_INT n, const cdouble a, cdouble *x) {
        cblas_zscal(n, &a, x, 1);
    }
    inline void scale(const MKL_INT n, const double a, cdouble *x) {
        cblas_zdscal(n, a, x, 1);
    }
} // MKL

/*
    ***********
    * Algebra *
    ***********
*/

template <class T, class U>
inline void init(T *vec, U size, T val){
    #pragma omp parallel for
    for (U i = 0; i < size; ++i) vec[i] = val;
}

template <class T>
inline void randInit(std::vector<T>& vec){
    std::mt19937 generatorD;
    std::uniform_real_distribution<double> distributionD(0.0,1.0);
    auto diceD = std::bind(distributionD,generatorD);
    for (auto& val : vec) val = diceD();
}

template <class T>
inline void copy(idx_t size, const T* vecIn, T* vecOut ) {
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) vecOut[i] = vecIn[i];
}

template <class T>
inline void copy(T* vecOut, T a, const T* vecIn, idx_t size) {
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) vecOut[i] = a * vecIn[i];
}

template <class T>
inline void axpy(T* y, T a, const T* x, idx_t size) {
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) y[i] += a * x[i];
}

template <class T>
inline void combine(T* y, T a1, const T* x1, T a2, const T* x2, idx_t size) {
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) y[i] = a1 * x1[i] + a2 * x2[i];
}

template <class T>
inline void scale(T* x, T a, idx_t size) {
    #pragma omp parallel for
    for (idx_t i = 0; i < size; ++i) x[i] *= a;
}

inline double dot(const VecD &v1, const VecD &v2) {
    assert_msg(v1.size()==v2.size(),"utils::dot, v1.size() != v2.size().");
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for(size_t i = 0; i < v1.size(); ++i){
        result += v1[i]*v2[i];
    }
    return result;
}

inline void normalize(VecD& v) {
    scale(v.data(),1.0/std::sqrt(dot(v,v)),v.size());
}

// vector dot prodoct
template <class T>
inline T mpiDot(const T* vconj, const T* v, idx_t size) {
    T partResult{0.0}, result{0.0};
    if constexpr (std::is_same<T, cdouble>::value) {
        double partResultReal = 0.0, partResultImag = 0.0;
        #pragma omp parallel for reduction(+:partResultReal,partResultImag)
        for (idx_t i = 0; i < size; ++i) {
            T tmp = std::conj(vconj[i]) * v[i];
            partResultReal += std::real(tmp);
            partResultImag += std::imag(tmp);
        }
        partResult = partResultReal + CPLX_I * partResultImag;
    } else if constexpr (std::is_same<T, double>::value) {
        #pragma omp parallel for reduction(+:partResult)
        for (idx_t i = 0; i < size; ++i) {
            partResult += vconj[i] * v[i];
        }
    }
    MPI_Allreduce(&partResult, &result, 1);
    return result;
}

template <typename T>
inline double mpiNorm(const T* v, idx_t size) {
    return std::sqrt(std::real(mpiDot(v, v, size)));
}
/*
    Eigen solver for a self adjoint matrix
*/
#endif
