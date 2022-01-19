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
#include "numbers"

using Basis_t = SpinBasis;

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
//    RedirectStdOut to("out/" + tostr(workerID) + ".out");

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
    int lx = *para.template get<int>("lx");
    int ly = *para.template get<int>("ly");
    auto nSite = (ly < 0) ? lx : lx * ly;
    int nu = *para.template get<int>("nu");
    int nd = *para.template get<int>("nd");
    double J1 = *para.template get<double>("J1");

    int specKD = *measurePara.template get<int>("spectraKrylovDim");
    int timeKD = *measurePara.template get<int>("timeKrylovDim");
    std::cout << std::setprecision(12);
    
    TriAngLattice latt(lx, ly, true);
    latt.addOrb({ORBITAL::SINGLE, 0, {0.0, 0.0, 0.0}}).construct();
    latt.print();
    for (int k = 0; k < nSite; ++k) {
        std::cout << "k = " << k << ", kxy = " << latt.getKxy(k) << std::endl;
    }

    Basis<Basis_t> b(&latt, nu, nd, 0, kidx, pidx);
    b.construct();
    std::cout << "nu = " << nu << ", nd = " << nd << std::endl;
    std::cout << "nu = " << BasisStateInterface::getNu() << ", nd = " << BasisStateInterface::getNd() << ", nPh = " << BasisStateInterface::getMaxPhPerSite() << std::endl;
    b.print();
    Link<dataType> J1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE},
                         J1, {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, -1.0, 1.0}});

    auto constructHam = [&](Hamiltonian<dataType, Basis_t> &h) {
        h.pushLinks({J1Link}).transform().construct();
    };
    timer.tik();
    Hamiltonian<dataType, Basis_t> H(&latt, &b, &b, true, true, 1, 1);
    constructHam(H);
    timer.print("H construction");
    H.printTrInteractions();
    H.printOnSiteInteraction();
    H.print("Hamiltonian info");
    // if (H.checkHermicity(0, 0)) {
    //     std::cout << "H is Hermitian!" << std::endl;
    // } else {
    //     std::cout << "H is not Hermitian!" << std::endl;
    // }
    auto nev = *measurePara.template get<int>("nev");
    H.diag(nev);

    double time = 0.0;
    VecD timePoints{time};

    if (opt(measurePara, "eqSpec")) {
        for (int kDyn = 0; kDyn < nSite; ++kDyn) {
            // int kDyn = *measurePara.template get<int>("kDynamic");
            int kf = (kidx == -1) ? -1 : (kidx + kDyn) % nSite;
            {
                Basis<Basis_t> bf(&latt, nu, nd, 0, kf, pidx);
                bf.construct();
                Hamiltonian<dataType, Basis_t> Hf(&latt, &bf, &bf, true, true, 1, 1);
                constructHam(Hf);

                if (opt(measurePara, "Skw")) {
                    OpK<Sz, Basis_t> szk(kDyn, ORBITAL::SINGLE, &latt, &b, &bf);
                    szk.construct();
                    SPECTRASolver<dataType> spectra(&Hf, H.getEval(), &szk, H.getEvec(), specKD);
                    spectra.compute();
                    if (isMaster) spectra.save(path.SkwDir + "/Lanczos/k" + tostr(kDyn), 0);
                    double skStatic = measureStaticStrucFact(&szk, H.getEvec());
                    if (isMaster) save(&skStatic, 1, path.SkwDir + "/k" + tostr(kDyn));
                    double dt2 = *measurePara.template get<double>("dt2");
                    int steps2 = *measurePara.template get<int>("steps2");
                    doubleTimeCorrelator(&szk, &H, &Hf, H.getEvec(), timeKD, dt2, steps2, path.SkwDir + "/time/k" + tostr(kDyn), isMaster);
                }
            }
        }
    }
    
    MPI_Finalize();
    return 0;
}
