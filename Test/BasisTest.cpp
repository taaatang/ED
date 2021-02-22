//
//  BasisTest.cpp
//  ED
//
//  Created by tatang on 8/4/19.
//  Copyright © 2019 tatang. All rights reserved.
//
#define OMP_
// #undef OMP_

#include "globalPara.hpp"
#include "PARPACKSolver.hpp"
#include "LANCZOSIterator.hpp"

#include <iostream>
#include <iomanip> // std::setprecision
#include <fstream>
#include <stdlib.h> // system
#include <chrono>
#include <mpi.h>
#ifdef OMP_
    #include <omp.h>
#endif

int main(int argc, const char * argv[]) {
    
}