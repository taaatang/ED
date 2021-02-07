#ifndef __GLOBALPARA_H__
#define __GLOBALPARA_H__

//
// globalPara.hpp
// ED
//
// Created by tatang on 11/19/19
//
#include "Global/globalType.hpp"

const bool DEBUG = false;
/*
    *******
    * MPI *
    *******
*/
const int MPI_MASTER = 0;

/*
    *****************
    * File I/O Path *
    *****************
*/
// cori
// const std::string ROOT_DATA_PATH = "/global/project/projectdirs/m2757/tatang";  
// sherlock
const std::string ROOT_DATA_PATH = "/oak/stanford/groups/tpd/tatang"; 
// acer
// const std::string ROOT_DATA_PATH = "/home/tatang/Project/DATA"; 

// const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/Photodoping";
// const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/TriAngtJ";
const std::string PROJECT_DATA_PATH = ROOT_DATA_PATH + "/TriAngHeis";
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
    * Geometry Class *
    ******************
*/
const int MAX_BOND_NUM = 6;
/*
    ******************
    * Spin Operators *
    ******************
*/
const int SPIN_DIM = 2;
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
/*
    *****************
    * PARPACKSolver *
    *****************
*/
// MAX parpack iterations
const int PARPACK_MAXITERATION = 5000;
//MIN and MAX ncv
const int PARPACK_MINNCV = 7;
// const int PARPACK_MAXNCV = 500;

/*
    *******************
    * LANCZOSIterator *
    *******************
*/
const LANCZOS_OPTION LANCZOS_DEFAULT_OPTION = ALPHA_BETA;
#endif // __GLOBALPARA_H__