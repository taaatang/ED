#include <iostream>
#include <sstream>
#include "Utils/paras.hpp"
#include "Global/config.hpp"
#include "Solver/PARPACKSolver.hpp"
#include "Solver/Spectra.hpp"
#include "Pulse/pulse.hpp"
#include "Solver/TimeEvolver.hpp"

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
        pulsePara.print(path.pumpFile);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ground state
    setBasics(para, latt, Bi, H);
    timer.tik();
    Bi->construct(opt(para, "basis"), path.getBasisDir(Bi->getkIndex(),  Bi->getPGIndex()));
    if (isMaster) Bi->print();
    H->construct();
    timer.tok();
    if (isMaster) {
        timer.print("Hamiltonian construction");
        // H->printLinks();
        H->print("Hamiltonian from master worker");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // cout<<setprecision(10);
    H->diag();
    cdouble w0 = H->getEval();
    cdouble* gstate = H->getEvec();

    if (measure("conductivity")) {
        // conductivity
        timer.tik();
        Current J(latt.get(), Bi.get());
        J.pushLinks(HubbardLink(*latt));
        auto ds = std::vector<std::string> {"x", "y", "z"};
        for (auto d : ds) {
            J.setDirection(d);
            // J.print();
            J.construct();
            SPECTRASolver<cdouble> spectra(H.get(), w0, &J, gstate, 400);
            spectra.compute();
            if (isMaster) {
                spectra.save(path.sigmaDir + "/" + d);
            }
        }
        timer.tok();
        if (isMaster) timer.print("Conductivity");
    }

    if (measure("Akw")) {
        // Akw
        timer.tik();
        auto occi = Bi->getOcc();
        auto ki = Bi->getkIndex();
        for (int kf = 0; kf < latt->getSiteNum(); ++kf) {
            for (auto spin : std::vector<SPIN> {SPIN::UP, SPIN::DOWN}) {
                for (auto pm : std::vector<LADDER> {LADDER::PLUS, LADDER::MINUS}) {
                    auto occf = occi;
                    auto& n = (spin == SPIN::UP) ? occf.at(0) : occf.at(1);
                    auto dn = (pm == LADDER::PLUS) ? 1 : -1;
                    n += dn;
                    setbasis(para, Bf, latt.get(), occf.at(0), occf.at(1), kf, -1);
                    setham(para, Hf, latt.get(), Bf.get());
                    Bf->construct(opt(para, "basis"), path.getBasisDir(Bf->getkIndex(),  Bf->getPGIndex()));
                    Hf->construct();
                    for (auto& orb : latt->getUnitCell()) {
                        CkOp<dataType> ck(latt.get(), Bi.get(), Bf.get());
                        ck.set(pm, spin, orb);
                        ck.construct();
                        SPECTRASolver<dataType> spectra(Hf.get(), w0, &ck, gstate, 400);
                        spectra.compute();
                        if (isMaster) {
                            std::ostringstream os;
                            os<<"/ki"<<tostr(ki)<<"/kf"<<tostr(kf)<<"/orb"<<orb.orbid<<"/"<<spin<<"_"<<pm;
                            spectra.save(path.AkwDir + os.str());
                        }
                    }
                }
            }
        }
        timer.tok();
        if (isMaster) timer.print("Single Particle Spectra");
    }

    if (measure("pump")) {
        // time evolution
        timer.tik();
        Pulse pulse;
        setpulse(pulsePara, pulse);
        if (isMaster) pulse.print();
        H->setPeierls(&pulse);
        // H->printPeierls();
        H->construct();
        timer.tok();
        if (isMaster) {
            timer.print("Set H(t)");
            H->print("H(t) from worker " + tostr(workerID));
        }

        int krylovDim = 15;
        TimeEvolver<cdouble> Tevol(gstate, H.get(), krylovDim);

        Nocc occ(latt.get(), Bi.get()); 
        occ.construct();

        timer.tik();
        if (isMaster) pulse.progressBar(100);
        // if (isMaster) std::cout<<"Progress  :";
        while (H->next()) {
            Tevol.evolve(pulse.getdt());
            occ.count(Tevol.getVec());
            if (isMaster) {
                pulse.progress();
                // if (pulse.getCount() % (pulse.getStepNum() / 100) == 0) std::cout<<">"<<std::flush;
            }
        }
        if (isMaster) std::cout<<"\nFinished\n";
        timer.tok();
        if (isMaster) timer.print("Timer evolution");

        if (isMaster) {
            occ.save(path.pumpDir);
        }
    }

    MPI_Finalize();
    return 0;
}