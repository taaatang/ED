//
//  Spectra.hpp
//  ED
//
//  Created by tatang on 2/20/20.
//  Copyright © 2019 tatang. All rights reserved.
//
#ifndef Spectra_hpp
#define Spectra_hpp

#include "LANCZOSIterator.hpp"

template <class T>
class SPECTRASolver: public LANCZOSIterator<T>{
private:
    // pointer to Hamiltonian and Operator A
    SparseMatrix<T> *A;
    // ground state energy
    cdouble w0;
    // initial state
    std::vector<T> vec;
public:
    SPECTRASolver(SparseMatrix<T> *H_, cdouble w0_, SparseMatrix<T> *A_, T* vec_, int krylovdim_):LANCZOSIterator<T>(H_,krylovdim_), \
    A(A_), w0(w0_) {
        assert_msg(A_->getDim() == H_->getColDim(), "Matrix dimension mismatch in SPECTRASolver!");
        vec.resize(A->getnlocmax());
        A->MxV(vec_,vec.data());
    }
    ~SPECTRASolver(){}
    void compute(){this->run(vec.data());}
    void save(std::string dataPath, int stateID = 0){
        mkdir_fs(dataPath);
        std::ofstream outfile;
        auto stateLabel =tostr(stateID);
        ::save<cdouble>(&w0, 1, &outfile, dataPath + "/w0_" + stateLabel);
        ::save<double>(&(LANCZOSIterator<T>::vecNorm), 1, &outfile, dataPath + "/vecNorm_" + stateLabel);
        ::save<double>((LANCZOSIterator<T>::alpha).data(), (int)((LANCZOSIterator<T>::alpha).size()), &outfile, dataPath + "/alpha_" + stateLabel);
        ::save<double>((LANCZOSIterator<T>::beta).data(), (int)((LANCZOSIterator<T>::beta).size()), &outfile, dataPath + "/beta_" + stateLabel);
        ::save<double>((LANCZOSIterator<T>::othErr).data(), (int)((LANCZOSIterator<T>::othErr).size()), &outfile, dataPath + "/othErr_" + stateLabel);
    };
};

#endif //Spectra_hpp

