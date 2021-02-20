#include <iostream>
#include <sstream>
#include "Utils/utils.hpp"
#include "Utils/paras.hpp"
#include "Global/config.hpp"
#include "Solver/PARPACKSolver.hpp"
#include "Solver/Spectra.hpp"
#include "Solver/TimeEvolver.hpp"

using namespace std;

int main( ) {

    MPI_Init(NULL, NULL);
    int workerID, workerNum;
    mpi_info(workerID, workerNum);
    Timer timer;

    auto isMaster = [=] {
        return (workerID == MPI_MASTER);
    };

    if (isMaster()) {
        path.make();
        modelPara.print(path.parameterFile);
        pulsePara.print(path.pumpFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ground state
    setBasics(modelPara, latt, Bi, H, workerID);
    // cout<<setprecision(10);
    H->solver.diag();
    cdouble* w0 = H->solver.getEigval();
    cdouble* gstate = H->solver.getEigvec();

    // conductivity
    Current J(latt.get(), Bi.get());
    J.pushLinks(HubbardLink());
    auto ds = std::vector<std::string> {"x", "y", "z"};
    for (auto d : ds) {
        J.setDirection(d);
        J.construct();
        SPECTRASolver<cdouble> spectra(H.get(), w0[0], &J, gstate, 400);
        spectra.compute();
        if (isMaster()) {
            spectra.saveData(path.sigmaDir + "/" + d);
        }
    }

    // Akw
    // auto occ = Bi->getOcc();
    // auto ki = Bi->getkIndex();
    // for (int kf = 0; kf < latt->getSiteNum(); ++kf) {
    //     for (auto spin : std::vector<SPIN> {SPIN::UP, SPIN::DOWN}) {
    //         for (auto pm : std::vector<LADDER> {LADDER::PLUS, LADDER::MINUS}) {
    //             auto occf = occ;
    //             auto& n = (spin == SPIN::UP) ? occf.at(0) : occf.at(1);
    //             auto dn = (pm == LADDER::PLUS) ? 1 : -1;
    //             n += dn;
    //             setbasis(modelPara, Bf, latt.get(), occf.at(0), occf.at(1), kf, -1);
    //             setham(modelPara, Hf, latt.get(), Bf.get());
    //             Bf->gen();
    //             Hf->construct();
    //             for (auto& orb : latt->getUnitCell()) {
    //                 CkOp<dataType> ck(latt.get(), Bi.get(), Bf.get());
    //                 ck.set(pm, spin, orb);
    //                 ck.construct();
    //                 SPECTRASolver<dataType> spectra(Hf.get(), w0[0], &ck, gstate, 400);
    //                 spectra.compute();
    //                 if (isMaster()) {
    //                     std::ostringstream os;
    //                     os<<"/ki"<<tostr(ki)<<"/kf"<<tostr(kf)<<"/"<<orb.orb<<"_"<<spin<<"_"<<pm;
    //                     spectra.saveData(path.AkwDir + os.str());
    //                 }
    //             }
    //         }
    //     }
    // }

    // // time evolution
    // timer.tik();
    // setpulse(pulsePara, pulse);
    // H->setPeierls(&pulse);
    // H->printPeierls();
    // H->construct();
    // timer.tok();
    // if (isMaster()) {
    //     timer.print("Set H(t)");
    //     H->print("H(t) from worker " + tostr(workerID));
    // }

    // int krylovDim = 15;
    // TimeEvolver<cdouble> Tevol(gstate, H.get(), krylovDim);

    // Nocc occ(latt.get(), Bi.get()); 
    // occ.construct();

    // timer.tik();
    // while (H->next()) {
    //     Tevol.evolve(pulse.getdt());
    //     occ.count(Tevol.getVec());
    // }
    // timer.tok();
    // if (isMaster()) timer.print("Timer evolution");

    // if (isMaster()) {
    //     occ.save(path.pumpDir);
    // }

    MPI_Finalize();
    return 0;
}