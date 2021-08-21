#pragma once

//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
#include "global/globalType.hpp"

/*
    *******
    * MPI *
    *******
*/
inline constexpr int MPI_MASTER = 0;

/*
    ******************
    * Math Constants *
    ******************
*/

inline constexpr double INFINITESIMAL = 1e-8;
inline constexpr double PI = 3.14159265358979323846;
inline constexpr std::complex<double> CPLX_I(0.0, 1.0);

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