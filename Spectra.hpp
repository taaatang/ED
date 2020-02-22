//
//  Spectra.hpp
//  ED
//
//  Created by tatang on 2/20/20.
//  Copyright © 2019 tatang. All rights reserved.
//
#ifndef Spectra_hpp
#define Spectra_hpp

#include "globalPara.hpp"
#include "Geometry.hpp"
#include "Basis.hpp"
#include "SparseMatrix.hpp"
#include "Operators.hpp"
#include "utils.hpp"
#include "LANCZOSIterator.hpp"
#include "HelperClass.hpp"

#include <iostream>
#include <fstream>
#include <stdlib.h> // system
#include <chrono>
#include <mpi.h>

template <class T>
class SPECTRASolver: public: LANCZOSIterator<T>{
private:
    // pointer to Hamiltonian and Operator A
    BaseMatrix<T> *H, *A;
    // ground state energy
    cdouble w0;
    // initial state
    std::vector<T> vec;
public:
    SPECTRASolver(BaseMatrix<T> *H_, double w0_, BaseMatrix<T> *A_, T* vec_, ind_int vecSize_, int krylovdim_):H(H_),A(A_),w0(w0_),vec(vec_),\
        LANCZOSIterator<T>(H_,krylovdim_){
        A->setBuf(vecSize_);
        vec.resize(H->get_nlocmax());
        A->MxV(vec_,vec.data());
        H->setBuf();
    }
    ~SPECTRASolver(){}
    void compute(){run(vec.data());}
    void save(std::string dataPath){
        std::ofstream outfile;
        save<cdouble>(&w0, 1, &outfile, dataPath + "/w0");
        save<cdouble>(&(LANCZOSIterator::vecNorm), 1, &outfile, dataPath + "/vecNorm");
        save<double>(LANCZOSIterator::alpha.data(), (int)lanczos.alpha.size(), &outfile, dataPath + "/alpha");
        save<double>(LANCZOSIterator::beta.data(), (int)lanczos.beta.size(), &outfile, dataPath + "/beta");
    };
};

#endif

