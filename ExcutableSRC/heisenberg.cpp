#include <iostream>
#include <sstream>
#include "Utils/paras.hpp"
#include "Global/config.hpp"
#include "Solver/PARPACKSolver.hpp"
#include "Solver/Spectra.hpp"

using namespace std;

int main( ) {

    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
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
    if (measure("basis") && !para.geti("basis")) {
        setlatt(para, latt);
        for (int kidx = workerID; kidx < latt->getSiteNum(); kidx += workerNum) {
            Basis b(para.getmodel(), latt.get(), {para.geti("nu"), para.geti("nd")}, kidx);
            b.construct();
            b.save(path.getBasisDir(kidx));
        }
    }

    // ground state
    cdouble* w0 = nullptr;
    cdouble* gstate = nullptr;
    if (measure("ground state")) {
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
        H->diag(5);
        w0 = H->getEval();
        gstate = H->getEvec();
    }
    // spin spin correlation
    if (measure("SiSj") && gstate) {
        timer.tik();
        SSOp<dataType> SS(latt.get(), Bi.get());
        std::vector<dataType> ssvals;
        for (int i = 0; i < latt->getSiteNum(); ++i) {
            SS.setr(i);
            SS.construct();
            auto val = SS.vMv(gstate, gstate) / dataType(latt->getSiteNum());
            if (isMaster) std::cout<<"ss(r="<<i<<")="<<val<<std::endl;
            ssvals.push_back(val);
        }
        if (isMaster) {
            save<dataType>(ssvals.data(), int(ssvals.size()), path.SiSjFile);
        }
        timer.print("Spin Spin Correlation");
    }

    // dynamic spin structure factor
    if (measure("Skw") && w0 && gstate) {
        timer.tik();
        auto occi = Bi->getOcc();
        for (int kf = 0; kf < latt->getSiteNum(); ++kf) {
            setbasis(para, Bf, latt.get(), occi.at(0), occi.at(1), kf, -1);
            setham(para, Hf, latt.get(), Bf.get());
            Bf->construct(opt(para,"basis"), path.getBasisDir(Bf->getkIndex(), Bf->getPGIndex()));
            Hf->construct();
            SzkOp<dataType> sk(latt.get(), Bi.get(), Bf.get());
            sk.construct();
            SPECTRASolver<dataType> spectra(Hf.get(), w0[0], &sk, gstate, 400);
            spectra.compute();
            if (isMaster) {
                spectra.saveData(path.SkwDir + "/kf" + tostr(kf));
            }
        }
        timer.print("Skw");
    }

    // Raman spectra
    if (measure("Raman") && w0 && gstate) {
        timer.tik();
        auto J1 = para.getd("J1");
        auto J2 = para.getd("J2");
        for (auto channel : std::vector<std::string> {"A1", "A2", "E21", "E22"}) {
            Hamiltonian<LATTICE_MODEL::HEISENBERG, dataType> R(latt.get(), Bi.get(), Bi.get());
            R.pushLinks(RamanChannel(channel, J1, J2, *latt));
            R.construct();
            SPECTRASolver<dataType> spectra(H.get(), w0[0], &R, gstate, 400);
            spectra.compute();
            if (isMaster) {
                spectra.saveData(path.RamanDir + "/" + channel);
            }
        }
        timer.print("Raman Spectra");
    }
    MPI_Finalize();
    return 0;
}