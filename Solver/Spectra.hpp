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
    SparseMatrix<T> *H, *A;
    // ground state energy
    cdouble w0;
    // initial state
    std::vector<T> vec;
public:
    SPECTRASolver(SparseMatrix<T> *H_, cdouble w0_, SparseMatrix<T> *A_, T* vec_, ind_int vecSize_, int krylovdim_):H(H_),A(A_),w0(w0_),\
        LANCZOSIterator<T>(H_,krylovdim_){
        std::cout<<"in spec init\n";
        A->setBuf(vecSize_);
        vec.resize(H->get_nlocmax());
        std::cout<<"begin A->MxV\n";
        A->MxV(vec_,vec.data());
        std::cout<<"finished MxV\n";
        H->setBuf();
        std::cout<<"finished spec init\n";
    }
    ~SPECTRASolver(){}
    void compute(){run(vec.data());}
    void saveData(std::string dataPath){
        std::ofstream outfile;
        save<cdouble>(&w0, 1, &outfile, dataPath + "/w0");
        save<cdouble>(&(LANCZOSIterator<T>::vecNorm), 1, &outfile, dataPath + "/vecNorm");
        save<double>((LANCZOSIterator<T>::alpha).data(), (int)((LANCZOSIterator<T>::alpha).size()), &outfile, dataPath + "/alpha");
        save<double>((LANCZOSIterator<T>::beta).data(), (int)((LANCZOSIterator<T>::beta).size()), &outfile, dataPath + "/beta");
    };
};

#endif //Spectra_hpp

