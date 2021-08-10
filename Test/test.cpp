//
// Created by tatang on 8/2/21.
//

#include "Basis/BasisState.hpp"
#include "Basis/Basis.hpp"
#include "Operator/localOperators.hpp"
#include "Operator/Operators.hpp"
#include "Solver/Spectra.hpp"
#include "Measure/config.hpp"
#include "Measure/measure.hpp"

using Basis_t = ElectronPhononBasis;

int main(int argc, char *argv[]) {
    std::string dataDir = "/Volumes/Sandisk/ProjectData/dmrg/benchmark/ED";
    int workerID{0}, workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(isMaster);

    assert_msg(argc > 1, "No inputDir!");
    std::string inputDir(argv[1]);
    Parameters pathPara(inputDir, {"path.txt"});
    Parameters para(inputDir, {"lattice.txt", "hamiltonian.txt"});
    Parameters pulsePara(inputDir, {"pulse.txt"});
    Parameters measurePara (inputDir, {"measure.txt"});
    Path path(&pathPara, &para, &pulsePara);
    if (isMaster) {
        path.make(measurePara);
        para.print(path.parameterFile);
    }

    ElectronBasis::setAllowDoubleOcc(true);
    int nSite = *para.template get<int>("lx");
    int nu = *para.template get<int>("nu");
    int nd = *para.template get<int>("nd");
    int maxPhPerSite = *para.template get<int>("maxPhPerSite");
    double tdp = *para.template get<double>("tdp");
    double Vd = *para.template get<double>("Vd");
    double Vp = *para.template get<double>("Vp");
    double Udd = *para.template get<double>("Udd"); 
    double Upp = *para.template get<double>("Upp");
    double gd = *para.template get<double>("gd");
    double gp = *para.template get<double>("gp");
    int specKD = *measurePara.template get<int>("spectraKrylovDim");
    int timeKD = *measurePara.template get<int>("timeKrylovDim");
    std::cout << std::setprecision(4);
    
    SquareLattice latt(nSite, 1, true);
    latt.addOrb({ORBITAL::Dx2y2, 0, {0.0, 0.0, 0.0}}).addOrb({ORBITAL::Px, 1, {0.5, 0.0, 0.0}}).construct();
    // latt.print();

    Basis<Basis_t> b(&latt, nu, nd, maxPhPerSite, -1);
    b.construct();
    // std::cout << b << std::endl;
    
    Hamiltonian<dataType, Basis_t> H(&latt, &b, &b, true, true, 1, 1);
    Link<dataType> tCuPx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px},    -tdp, {{0.5, 0.0, 0.0}});
    Link<dataType> tPxCu(LINK_TYPE::HOPPING_T, {ORBITAL::Px,    ORBITAL::Dx2y2},  tdp, {{0.5, 0.0, 0.0}});
    Link<dataType> gCu(LINK_TYPE::NCHARGE_SITE_PHONON, {ORBITAL::Dx2y2, ORBITAL::Dx2y2}, gd, {{0.0, 0.0, 0.0}});
    Link<dataType> gPx(LINK_TYPE::NCHARGE_SITE_PHONON, {ORBITAL::Px, ORBITAL::Px}, gp, {{0.0, 0.0, 0.0}});
    H.pushLinks({tCuPx, tPxCu, gCu, gPx}).pushV({ORBITAL::Dx2y2}, Vd).pushV({ORBITAL::Px}, Vp).pushU({ORBITAL::Dx2y2}, Udd).pushU({ORBITAL::Px}, Upp).transform().construct();
    // H.printLinks();
    H.diag();
    // H.print("Hamiltonian info", std::cout, false);
//    std::cout << "eval: " << H.getEval() << std::endl;

    NelOp<Basis_t> nel(&latt, &b, &b);
    nel.construct();
    nel.count(H.getEvec());
    if (isMaster) std::cout << "charge distribution: " << nel.lastCount() << std::endl;

    NphOp<Basis_t> nph(&latt, &b, &b);
    nph.construct();
    nph.count(H.getEvec());
    if (isMaster) std::cout << "phonon distribution: " << nph.lastCount() << std::endl;

    double time = 0.0;
    VecD timePoints{time};

    if (opt(measurePara, "conductivity")) {
        Current<Basis_t> jx(&latt, &b, &b);
        jx.pushLinks({tCuPx, tPxCu});
        jx.setDirection("x");
        jx.construct();
        // jx.print("x current");
        SPECTRASolver<cdouble> spectra(&H, H.getEval(), &jx, H.getEvec(), specKD);
        spectra.compute();
        if (isMaster) {
            spectra.save(path.sigmaDir + "/x", 0);
        }
    }

    ORBITAL orb = (*measurePara.template get<int>("skOrbId") == 0) ? ORBITAL::Dx2y2 : ORBITAL::Px;
    std::string orbName = (orb == ORBITAL::Dx2y2) ? "Cu" : "Ox";
    std::vector<std::unique_ptr<SzkOp<dataType, Basis_t>>> szk;
    for (int i = 0; i < nSite; ++i) {
        szk.push_back(std::unique_ptr<SzkOp<dataType, Basis_t>>(new SzkOp<dataType, Basis_t>(i, orb, &latt, &b, &b)));
        szk[i]->construct();
    }
    // std::vector<dataType> szkVal;
    // for (int i = 0; i < nSite; ++i) {
    //     std::vector<dataType> v1(szk[i]->getnloc(), 0.0);
    //     szk[i]->MxV(H.getEvec(), v1.data());
    //     szkVal.push_back(mpiDot(v1.data(), v1.data(), v1.size()));
    //     // szkVal.push_back(szk[(nSite - i) % nSite]->vMv(H.getEvec(), v1.data()));
    // }
    // std::cout << "szk:\n" << szkVal << std::endl;
    if (opt(measurePara, "Skw")) {
        for (int k = 0; k < szk.size(); ++k) {
            SPECTRASolver<dataType> spectra(&H, H.getEval(), szk[k].get(), H.getEvec(), 100);
            spectra.compute();
            if (isMaster) spectra.save(path.SkwDir + "/" + orbName + "/Lanczos/k" + tostr(k), 0);
            doubleTimeCorrelator(szk[(nSite - k) % nSite].get(), szk[k].get(), &H, H.getEvec(), timeKD, 0.5, 800, path.SkwDir + "/" + orbName + "/time/k" + tostr(k), isMaster);
        }
    }
    
    if (opt(measurePara, "pump")) {
        Pulse pulse;
        setpulse(pulsePara, pulse);
        H.setPeierls(&pulse);
        H.construct();
        TimeEvolver<dataType> Tevol(H.getEvec(), &H, timeKD);
        if (isMaster) pulse.progressBar(20);
        timer.tik();
        while (H.next()) {
            Tevol.evolve(pulse.getdt());
            nel.count(Tevol.getVec());
            nph.count(Tevol.getVec());
            time += pulse.getdt();
            timePoints.push_back(time);
            if (isMaster) {
                pulse.progress();
            }
        }
        timer.print("pump time evolution");
        H.turnOffPulse();

        if (opt(measurePara, "neqSkw")) {
            double dt1 = *measurePara.template get<double>("dt1");
            int steps1 = *measurePara.template get<int>("steps1");
            double dt2 = *measurePara.template get<double>("dt2");
            int steps2 = *measurePara.template get<int>("steps2");
            for (int i = 0; i < steps1; ++i) {
                nel.count(Tevol.getVec());
                nph.count(Tevol.getVec());
                time += dt1;
                timePoints.push_back(time);
                int k = 2;
                // for (int k = 0; k < szk.size(); ++k) {
                timer.tik();
                doubleTimeCorrelator(szk[(nSite - k) % nSite].get(), szk[k].get(), &H, Tevol.getVec(), timeKD, dt2, steps2, path.NeqSkwDir + "/" + orbName + "/k" + tostr(k) + "/t1_" + tostr(i), isMaster);
                // }
                Tevol.evolve(dt1);
                timer.print("szk(t1,t2)");
            }
        }
        if (isMaster) {
            nel.save(path.NeqNelDir);
            nph.save(path.NeqNphDir);
            save<double>(timePoints.data(), int(timePoints.size()), path.NeqNelDir + "/timePoints");
        }
    }
    
    MPI_Finalize();
    return 0;
}
