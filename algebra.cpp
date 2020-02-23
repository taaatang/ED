//
// algebra.hpp
// ED
//
// Created by tatang on 1/15/2020
//
//

#include "algebra.hpp"

void diagTri(std::vector<double>* a, std::vector<double>* b, std::vector<double>* U){
    // std::cout<<"alpha size:"<<a->size()<<", beta size:"<<b->size()<<std::endl;
    int dim = a->size();
    assert(b->size()>=(dim-1));
    assert(U->size()>=(dim*dim));
    // double symmetric triangular eigen divide and conquer algorithm
    LAPACKE_dstedc(LAPACK_COL_MAJOR, 'I', dim, a->data(), b->data(), U->data(), dim);
}