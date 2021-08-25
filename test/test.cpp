//
// Created by tatang on 8/2/21.
//

#include "basis/basisState.hpp"
#include "basis/basis.hpp"
#include "operator/localOperators.hpp"
#include "operator/operators.hpp"
#include "solver/spectra.hpp"
#include "measure/config.hpp"
#include "measure/measure.hpp"
#include "utils/progressBar.hpp"

using Basis_t = ElectronPhononBasis;

int main(int argc, char *argv[]) {
    int workerID{0};
    int workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(true);
    if (isMaster) {
        mkdir_fs("out");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    RedirectStdOut to("out/" + tostr(workerID) + ".out");

    assert_msg(argc > 1, "No inputDir!");
    std::string inputDir(argv[1]);
    Parameters pathPara(inputDir, {"path.txt"});
    Parameters para(inputDir, {"lattice.txt", "hamiltonian.txt"});
    Parameters pulsePara(inputDir, {"pulse.txt"});
    Parameters measurePara(inputDir, {"measure.txt"});

    Path path(&pathPara, &para, &pulsePara);
    if (isMaster) {
        path.make(measurePara);
        para.print(path.parameterFile);
    }

    ElectronBasis::setAllowDoubleOcc(true);
    int kidx = *para.template get<int>("kidx");
    int pidx = *para.template get<int>("pidx");
    int nSite = *para.template get<int>("lx");
    int nu = *para.template get<int>("nu");
    int nd = *para.template get<int>("nd");
    int maxPhPerSite = *para.template get<int>("maxPhPerSite");
    double tdp = *para.template get<double>("tdp");
    double Vd = *para.template get<double>("Vd");
    double Vp = *para.template get<double>("Vp");
    double Udd = *para.template get<double>("Udd"); 
    double Upp = *para.template get<double>("Upp");
    double wd = *para.template get<double>("wd");
    double wp = *para.template get<double>("wp");
    double wdwp = std::sqrt(wd * wp) / 2.0;
    double gd = *para.template get<double>("gd");
    double gp = *para.template get<double>("gp");
    int specKD = *measurePara.template get<int>("spectraKrylovDim");
    int timeKD = *measurePara.template get<int>("timeKrylovDim");
    std::cout << std::setprecision(12);
    
    SquareLattice latt(nSite, 1, true);
    latt.addOrb({ORBITAL::Dx2y2, 0, {0.0, 0.0, 0.0}}).addOrb({ORBITAL::Px, 1, {0.5, 0.0, 0.0}}).construct();
    // latt.print();

    Basis<Basis_t> b(&latt, nu, nd, maxPhPerSite, kidx, pidx);
    b.construct();
    b.print();
    Hamiltonian<dataType, Basis_t> H(&latt, &b, &b, true, true, 1, 1);
    Link<dataType> tCuPx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px},    -tdp, {{0.5, 0.0, 0.0}});
    Link<dataType> tPxCu(LINK_TYPE::HOPPING_T, {ORBITAL::Px,    ORBITAL::Dx2y2},  tdp, {{0.5, 0.0, 0.0}});
    Link<dataType> kCuPx(LINK_TYPE::XiXj,      {ORBITAL::Dx2y2, ORBITAL::Px},    wdwp, {{0.5, 0.0, 0.0}});
    Link<dataType> kPxCu(LINK_TYPE::XiXj,      {ORBITAL::Px,    ORBITAL::Dx2y2}, wdwp, {{0.5, 0.0, 0.0}});
    Link<dataType> gCu(LINK_TYPE::NCHARGE_SITE_PHONON, {ORBITAL::Dx2y2, ORBITAL::Dx2y2}, gd, {{0.0, 0.0, 0.0}});
    Link<dataType> gPx(LINK_TYPE::NCHARGE_SITE_PHONON, {ORBITAL::Px, ORBITAL::Px}, gp, {{0.0, 0.0, 0.0}});
    timer.tik();
    H.pushLinks({tCuPx, tPxCu, gCu, gPx, kCuPx, kPxCu}).pushV({ORBITAL::Dx2y2}, Vd).pushV({ORBITAL::Px}, Vp).pushU({ORBITAL::Dx2y2}, Udd).pushU({ORBITAL::Px}, Upp).pushPhW0({ORBITAL::Dx2y2}, wd).pushPhW0({ORBITAL::Px}, wp).transform().construct();
    // H.pushLinks({tCuPx, tPxCu}).pushV({ORBITAL::Dx2y2}, Vd).pushV({ORBITAL::Px}, Vp).pushU({ORBITAL::Dx2y2}, Udd).pushU({ORBITAL::Px}, Upp).transform().construct();
    timer.print("H construction");
    H.printTrInteractions();
    H.print("Hamiltonian info");
    H.diag();

    NelOp<Basis_t> nel(&latt, &b, &b);
    nel.construct();
    nel.count(H.getEvec());
    std::cout << "charge distribution: " << nel.lastCount() << std::endl;

    NphOp<Basis_t> nph(&latt, &b, &b);
    nph.construct();
    nph.count(H.getEvec());
    std::cout << "phonon distribution: " << nph.lastCount() << std::endl;

    // for (auto orb : std::vector<ORBITAL>{ORBITAL::Dx2y2, ORBITAL::Px}) {
    //     OpK<NCharge<SPIN::UP>, Basis_t> n(0, orb, &latt, &b, &b);
    //     n.construct();
    //     auto val = n.vMv(H.getEvec(), H.getEvec());
    //     if (isMaster) std::cout << orb << " " << SPIN::UP << " has " << val << " charge" << std::endl;
    // }
    // for (auto orb : std::vector<ORBITAL>{ORBITAL::Dx2y2, ORBITAL::Px}) {
    //     OpK<NCharge<SPIN::DOWN>, Basis_t> n(0, orb, &latt, &b, &b);
    //     n.construct();
    //     auto val = n.vMv(H.getEvec(), H.getEvec());
    //     if (isMaster) std::cout << orb << " " << SPIN::DOWN << " has " << val << " charge" << std::endl;
    // }
    // for (auto orb : std::vector<ORBITAL>{ORBITAL::Dx2y2, ORBITAL::Px}) {
    //     OpK<NPhonon, Basis_t> n(0, orb, &latt, &b, &b);
    //     n.construct();
    //     auto val = n.vMv(H.getEvec(), H.getEvec());
    //     if (isMaster) std::cout << orb << " " << " has " << val << " phonon" << std::endl;
    // }
    double time = 0.0;
    VecD timePoints{time};

    if (opt(measurePara, "conductivity")) {
        Current<Basis_t> jx("x", &latt, &b, &b);
        jx.pushLinks({tCuPx, tPxCu});
        jx.setDirection();
        jx.transform();
        jx.construct();
        jx.print("x current");
        SPECTRASolver<cdouble> spectra(&H, H.getEval(), &jx, H.getEvec(), specKD);
        spectra.compute();
        if (isMaster) {
            spectra.save(path.sigmaDir + "/x", 0);
        }
    }

    for (int kDyn = 0; kDyn < nSite; ++kDyn) {
        ORBITAL orb = (*measurePara.template get<int>("skOrbId") == 0) ? ORBITAL::Dx2y2 : ORBITAL::Px;
        std::string orbName = (orb == ORBITAL::Dx2y2) ? "Cu" : "Ox";
        // int kDyn = *measurePara.template get<int>("kDynamic");
        int kf = (kidx == -1) ? -1 : (kidx + kDyn) % nSite;
        Basis<Basis_t> bf(&latt, nu, nd, maxPhPerSite, kf, pidx);
        bf.construct();
        Hamiltonian<dataType, Basis_t> Hf(&latt, &bf, &bf, true, true, 1, 1);
        Hf.pushLinks({tCuPx, tPxCu, gCu, gPx}).pushV({ORBITAL::Dx2y2}, Vd).pushV({ORBITAL::Px}, Vp).pushU({ORBITAL::Dx2y2}, Udd).pushU({ORBITAL::Px}, Upp).pushPhW0({ORBITAL::Dx2y2}, wd).pushPhW0({ORBITAL::Px}, wp).transform().construct();
        SzkOp<dataType, Basis_t> szk(kDyn, orb, &latt, &b, &bf);
        szk.transform();
        szk.construct();

        double dt2 = *measurePara.template get<double>("dt2");
        int steps2 = *measurePara.template get<int>("steps2");
        if (opt(measurePara, "Skw")) {
            SPECTRASolver<dataType> spectra(&Hf, H.getEval(), &szk, H.getEvec(), 100);
            spectra.compute();
            if (isMaster) spectra.save(path.SkwDir + "/" + orbName + "/Lanczos/k" + tostr(kDyn), 0);
            doubleTimeCorrelator(&szk, &H, &Hf, H.getEvec(), timeKD, dt2, steps2, path.SkwDir + "/" + orbName + "/time/k" + tostr(kDyn), isMaster);
        }
    }
    
    if (opt(measurePara, "pump")) {
        ORBITAL orb = (*measurePara.template get<int>("skOrbId") == 0) ? ORBITAL::Dx2y2 : ORBITAL::Px;
        std::string orbName = (orb == ORBITAL::Dx2y2) ? "Cu" : "Ox";
        int kDyn = *measurePara.template get<int>("kDynamic");
        int kf = (kidx == -1) ? -1 : (kidx + kDyn) % nSite;
        Basis<Basis_t> bf(&latt, nu, nd, maxPhPerSite, kf, pidx);
        bf.construct();
        Hamiltonian<dataType, Basis_t> Hf(&latt, &bf, &bf, true, true, 1, 1);
        std::cout << "construc Hf" << std::endl;
        Hf.pushLinks({tCuPx, tPxCu, gCu, gPx}).pushV({ORBITAL::Dx2y2}, Vd).pushV({ORBITAL::Px}, Vp).pushU({ORBITAL::Dx2y2}, Udd).pushU({ORBITAL::Px}, Upp).pushPhW0({ORBITAL::Dx2y2}, wd).pushPhW0({ORBITAL::Px}, wp).transform().construct();
        SzkOp<dataType, Basis_t> szk(kDyn, orb, &latt, &b, &bf);
        szk.transform();
        szk.construct();
        double dt2 = *measurePara.template get<double>("dt2");
        int steps2 = *measurePara.template get<int>("steps2");

        Pulse pulse;
        setpulse(pulsePara, pulse);
        pulse.profile();
        H.setPeierls(&pulse);
        H.transform();
        // H.printPeierls();
        // H.printTrInteractions();
        H.construct();

        // if (H.checkHermicity(1, 2)) {
        //     std::cout << "H is Hermitian!" << std::endl;
        // } else {
        //     std::cout << "H is not Hermitian!" << std::endl;
        // }
        H.print("H(t) info");

        TimeEvolver<dataType> Tevol(H.getEvec(), &H, timeKD);
        MPI_Barrier(MPI_COMM_WORLD);
        ProgressBar bar("pump", pulse.getStepNum(), 20, isMaster);
        timer.tik();
        VecD waveNorms;
        waveNorms.push_back(mpiNorm(Tevol.getVec(), H.getnloc()));
        while (H.next()) {
            Tevol.evolve(pulse.getdt());
            waveNorms.push_back(mpiNorm(Tevol.getVec(), H.getnloc()));
            nel.count(Tevol.getVec());
            nph.count(Tevol.getVec());
            time += pulse.getdt();
            timePoints.push_back(time);
            bar.progress();
        }
        timer.print("pump time evolution");
        H.turnOffPulse();
        double dt1 = *measurePara.template get<double>("dt1");
        int steps1 = *measurePara.template get<int>("steps1");
        for (int i = 0; i < steps1; ++i) {
            if (opt(measurePara, "neqOcc")) {
                waveNorms.push_back(mpiNorm(Tevol.getVec(), H.getnloc()));
                nel.count(Tevol.getVec());
                nph.count(Tevol.getVec());
            }
            if (opt(measurePara, "neqSkw")) {
                timer.tik();
                doubleTimeCorrelator(&szk, &H, &Hf, Tevol.getVec(), timeKD, dt2, steps2, path.NeqSkwDir + "/" + orbName + "/k" + tostr(kDyn) + "/t1_" + tostr(i), isMaster);
                timer.print("szk(t1,t2)");
                std::cout << std::endl;
            }
            time += dt1;
            timePoints.push_back(time);
            Tevol.evolve(dt1);
        }
    
        if (isMaster) {
            nel.save(path.NeqNelDir);
            nph.save(path.NeqNphDir);
            save<double>(timePoints.data(), timePoints.size(), path.NeqNelDir + "/timePoints");
            save<double>(waveNorms.data(), waveNorms.size(), path.NeqNelDir + "/waveNorms");
        }
    }
    
    MPI_Finalize();
    return 0;
}
