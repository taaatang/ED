#pragma once

//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
#include "Global/globalType.hpp"

constexpr bool OPT_USE_MKL = true;

constexpr bool DEBUG = false;
/*
    *******
    * MPI *
    *******
*/
constexpr int MPI_MASTER = 0;

/*
    ******************
    * Math Constants *
    ******************
*/

const double INFINITESIMAL = 1e-8;
const double PI = 3.14159265358979323846;
const std::complex<double> CPLX_I(0.0, 1.0);

/*
    ******************
    * Spin Operators *
    ******************
*/
constexpr int SPIN_DIM = 2;
/*
    **********************************
    * Sparse Matrix PARTITION Method *
    **********************************
*/
#ifdef DISTRIBUTED_STATE
    #define SPM_COL_PARTITION
    const MATRIX_PARTITION PARTITION = COL_PARTITION;
#else
    #define SPM_ROW_PARTITION
    const MATRIX_PARTITION PARTITION = ROW_PARTITION;
#endif