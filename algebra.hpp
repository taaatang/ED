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
#define _MKL_
#include "globalPara.hpp"
#define MKL_INT ind_int
#define MKL_Complex16 cdouble
#include "mkl.h"
#include "omp.h"
#include <vector>
/*
    LARPACK eigen solver for a real symmetric tridiagonal matrix
*/
namespace MKL{
    void diagTri(std::vector<double>* a, std::vector<double>* b, std::vector<double>* U);

    void MxV(const MKL_INT& rows, const MKL_INT& cols, const MKL_Complex16& alpha, const std::vector<MKL_Complex16>& val, const std::vector<MKL_INT>& colIndx, const std::vector<MKL_INT>& rowIndx , const MKL_Complex16 *vecIn , MKL_Complex16 *vecOut);

    void create(sparse_matrix_t& A, ind_int rows, ind_int cols, std::vector<MKL_INT>& rowInitList, std::vector<MKL_INT>& colList, std::vector<cdouble>& valList, MKL_INT mvNum=100);
    void destroy(sparse_matrix_t& A);
    void MxV2(sparse_matrix_t& A, MKL_Complex16* vin, MKL_Complex16* vout, MKL_Complex16 alpha=1.0, MKL_Complex16 beta=1.0);
}


/*
    Eigen solver for a self adjoint matrix
*/
#endif
