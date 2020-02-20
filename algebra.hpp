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
#include "globalPara.hpp"
#include "mkl.h"
#include "omp.h"
#include <vector>
/*
    LARPACK eigen solver for a real symmetric tridiagonal matrix
*/
void diagTri(std::vector<double>* a, std::vector<double>* b, std::vector<double>* U);

/*
    Eigen solver for a self adjoint matrix
*/
#endif
