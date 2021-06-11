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
#include "BiCGSTAB.hpp"

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

class SPECTRASolverBiCGSTAB {
public:
    SPECTRASolverBiCGSTAB(BaseMatrix<cdouble> *Ham, BaseMatrix<cdouble> *Op, cdouble *vec_, double w0, double wmin = 0.0, double wmax = 5.0, double dw = 0.01, double epsilon = 0.02);
    void compute();
    void save(std::string dataPath, int stateID = 0);
private:
    BaseMatrix<cdouble> *H;
    BaseMatrix<cdouble> *O;
    std::vector<cdouble> b;
    double w0;
    double wmin;
    double wmax;
    double dw;
    double epsilon;
    std::vector<double> spectra;
    std::vector<double> spectraInv;
};

SPECTRASolverBiCGSTAB::SPECTRASolverBiCGSTAB(BaseMatrix<cdouble> *Ham, BaseMatrix<cdouble> *Op, cdouble *vec_, double w0_, double wmin_, double wmax_, double dw_, double epsilon_) {
    H = Ham;
    O = Op;
    b.resize(O->getnlocmax());
    O->MxV(vec_, b.data());
    w0 = w0_;
    wmin = wmin_;
    wmax = wmax_;
    dw = dw_;
    epsilon = epsilon_;
}

void SPECTRASolverBiCGSTAB::compute() {
    spectra.clear();
    spectraInv.clear();
    Identity<cdouble> eye(H->getDim());
    for (auto w = wmin; w < wmax; w += dw) {
        cdouble z{-(w + w0), -epsilon};
        std::vector<cdouble> x(H->getnlocmax(), 0.0);
        int iterCount;
        double res;
        BiCGSTAB(H, z, b.data(), x.data(), H->getnloc(), &eye, iterCount, res);
        std::cout << "iter: " << iterCount << ", err: " << res << '\n';
        spectra.push_back(std::imag(mpiDot(b.data(), x.data(), H->getnloc())) / PI);
    }

    for (auto w = wmin; w < wmax; w += dw) {
        cdouble z{w - w0, -epsilon};
        std::vector<cdouble> x(H->getnlocmax(), 0.0);
        int iterCount;
        double res;
        BiCGSTAB(H, z, b.data(), x.data(), H->getnloc(), &eye, iterCount, res);
        std::cout << "iter: " << iterCount << ", err: " << res << '\n';
        spectraInv.push_back(std::imag(mpiDot(b.data(), x.data(), H->getnloc())) / PI);
    }
}

void SPECTRASolverBiCGSTAB::save(std::string dataPath, int stateID) {
    mkdir_fs(dataPath);
    std::ofstream outfile;
    auto stateLabel =tostr(stateID);
    ::save<double>(&w0, 1, &outfile, dataPath + "/w0_" + stateLabel);
    ::save<double>(&wmin, 1, &outfile, dataPath + "/wmin_" + stateLabel);
    ::save<double>(&wmax, 1, &outfile, dataPath + "/wmax_" + stateLabel);
    ::save<double>(&dw, 1, &outfile, dataPath + "/dw_" + stateLabel);
    ::save<double>(spectra.data(), int(spectra.size()), &outfile, dataPath + "/spectra_" + stateLabel);
    ::save<double>(spectraInv.data(), int(spectraInv.size()), &outfile, dataPath + "/spectraInv_" + stateLabel);
}

#endif //Spectra_hpp

