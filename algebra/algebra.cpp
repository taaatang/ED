//
// algebra.hpp
// ED
//
// Created by tatang on 1/15/2020
//
//

#include "algebra.hpp"

/*
    MKL diagonalize a tridiagonal matrix.
    input: 
        a: diagonal elemants
        b: subdiagonal elements
    output: 
        a: eigenvalues
        U(COL_MAJOR): eigen vectors
*/
void MKL::diagTri(std::vector<double>* a, std::vector<double>* b, std::vector<double>* U){
    // std::cout<<"alpha size:"<<a->size()<<", beta size:"<<b->size()<<std::endl;
    int dim = a->size();
    assert((int)b->size()>=(dim-1));
    assert((int)U->size()>=(dim*dim));
    // double symmetric triangular eigen divide and conquer algorithm
    LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', dim, a->data(), b->data(), U->data(), dim);
}

void MKL::create(sparse_matrix_t& A, idx_t rows, idx_t cols, std::vector<MKL_INT>& rowInitList, std::vector<MKL_INT>& colList, std::vector<MKL_Complex16>& valList, MKL_INT mvNum){
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_z_create_csr(&A, SPARSE_INDEX_BASE_ZERO, rows, cols, rowInitList.data(), rowInitList.data()+1, colList.data(), valList.data());
    mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descrA, mvNum);
    mkl_sparse_set_memory_hint(A, SPARSE_MEMORY_NONE);
    mkl_sparse_optimize(A);
}

void MKL::MxV(sparse_matrix_t& A, MKL_Complex16* vin, MKL_Complex16* vout, MKL_Complex16 alpha, MKL_Complex16 beta){
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descrA, vin, beta, vout);
}

void MKL::create(sparse_matrix_t& A, idx_t rows, idx_t cols, std::vector<MKL_INT>& rowInitList, std::vector<MKL_INT>& colList, std::vector<double>& valList, MKL_INT mvNum){
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, rows, cols, rowInitList.data(), rowInitList.data()+1, colList.data(), valList.data());
    mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descrA, mvNum);
    mkl_sparse_set_memory_hint(A, SPARSE_MEMORY_NONE);
    mkl_sparse_optimize(A);
}

void MKL::MxV(sparse_matrix_t& A, double* vin, double* vout, double alpha, double beta){
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descrA, vin, beta, vout);
}

void MKL::destroy(sparse_matrix_t& A){mkl_sparse_destroy(A);}

