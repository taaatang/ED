#include <iostream>
#include <sstream>
#include "utils/paras.hpp"
#include "global/config.hpp"
#include "solver/parpackSolver.hpp"
#include "solver/spectra.hpp"

int main( ) {

    MPI_Init(NULL, NULL);
    int workerID{0}; 
    int workerNum{0};
    mpiInfo(workerID, workerNum);
    bool isMaster = (workerID == MPI_MASTER);
    Timer timer(isMaster);

    auto measure = [&] (std::string key) {
        return opt(measurePara, key);
    };

    if (isMaster) {
        path.make();
        para.print(path.parameterFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // compute and save basis with translation symm
    if (measure("basis") && !*para.get<bool>("basis")) {
        setlatt(para, latt);
        for (int kidx = workerID; kidx < latt->getSiteNum(); kidx += workerNum) {
            Basis b(para.getmodel(), latt.get(), {*para.get<int>("nu"), *para.get<int>("nd")}, kidx);
            b.construct();
            b.save(path.getBasisDir(kidx));
        }
    }

    // eigenstate
    std::vector<dataType> w{};
    std::vector<dataType*> state{};
    int nstate = 0;
    if (measure("state")) {
        timer.tik();
        setBasics(para, latt, Bi, H);
        Bi->construct(opt(para, "basis"), path.getBasisDir(Bi->getkIndex(),  Bi->getPGIndex()));
        if (isMaster) Bi->print();
        H->construct();
        timer.print("Hamiltonian construction");
        if (isMaster) {
            H->printLinks();
            H->print("Hamiltonian from master worker");
        }
        // cout<<setprecision(10);
        auto nev = *measurePara.get<int>("nev");
        H->diag(nev);
        for (int n = 0; n < nev; ++n) {
            w.push_back(H->getEval(n));
            state.push_back(H->getEvec(n));
        }
        nstate = sort(w, state, measure("zeroTmp"));
        if (isMaster) {
            save(w.data(), (int)(w.size()), path.evalFile);
        }
        if (isMaster) {
            std::cout << "sorted evals:" << w << '\n';
        }
    }
    int krylovDim = *measurePara.get<int>("krylovDim");
    // spin spin correlation
    if (measure("SiSj")) {
        timer.tik();
        SSOp<dataType> SS(latt.get(), Bi.get());
        std::vector<std::vector<dataType>> ssvals;
        ssvals.resize(nstate);
        for (int i = 0; i < latt->getSiteNum(); ++i) {
            SS.setr(i);
            SS.construct();
            for (int n = 0; n < nstate; ++n) {
                auto val = SS.vMv(state[n], state[n]) / dataType(latt->getSiteNum());
                if (isMaster) std::cout<<"ss(r="<<i<<")="<<val<<std::endl;
                ssvals[n].push_back(val);
            }
        }
        for (int n = workerID; n < nstate; n += workerNum) {
            save(ssvals[n].data(), (int)ssvals[n].size(), path.SiSjFile + "_" + tostr(n));
        }
        timer.print("Spin Spin Correlation");
    }

    // dynamic spin structure factor
    if (measure("Skw")) {
        timer.tik();
        auto occi = Bi->getOcc();
        for (int kf = 0; kf < latt->getSiteNum(); ++kf) {
            setbasis(para, Bf, latt.get(), occi.at(0), occi.at(1), kf, -1);
            setham(para, Hf, latt.get(), Bf.get());
            Bf->construct(opt(para,"basis"), path.getBasisDir(Bf->getkIndex(), Bf->getPGIndex()));
            Hf->construct();
            SzkOp<dataType> sk(latt.get(), Bi.get(), Bf.get());
            sk.construct();
            for (int n = 0; n < nstate; ++n) {
                SPECTRASolver<dataType> spectra(Hf.get(), w[n], &sk, state[n], krylovDim);
                spectra.compute();
                if (isMaster) spectra.save(path.SkwDir + "/kf" + tostr(kf), n);
            }
        }
        timer.print("Skw");
    }

    // Raman spectra
    if (measure("Raman")) {
        timer.tik();
        auto J1 = *para.get<double>("J1");
        auto J2 = *para.get<double>("J2");
        for (auto channel : std::vector<std::string> {"A1", "A2", "E21", "E22"}) {
            Hamiltonian<MODEL::HEISENBERG, dataType> R(latt.get(), Bi.get(), Bi.get());
            R.pushLinks(RamanChannel(channel, J1, J2, *latt));
            R.construct();
            for (int n = 0; n < nstate; ++n) {
                SPECTRASolver<dataType> spectra(H.get(), w[n], &R, state[n], krylovDim);
                spectra.compute();
                if (isMaster) spectra.save(path.RamanDir + "/" + channel, n);
            }
        }
        timer.print("Raman Spectra");
    }
    MPI_Finalize();
    return 0;
}