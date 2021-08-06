//
// Created by tatang on 8/2/21.
//

#include "Basis/BasisState.hpp"
#include "Basis/Basis.hpp"
#include "Operator/localOperators.hpp"
#include "Operator/Operators.hpp"

using Basis_t = ElectronPhononBasis;

int main() {
    int workerID{0}, workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(isMaster);
    std::cout << std::setprecision(4);
    int nSite = 4;
    int nu = 2;
    int nd = 2;
    int nPho = 0;
    ElectronBasis::setAllowDoubleOcc(true);

    SquareLattice latt(nSite, 1, true);
    latt.addOrb({ORBITAL::Dx2y2, 0, {0.0, 0.0, 0.0}}).addOrb({ORBITAL::Px, 1, {0.5, 0.0, 0.0}});
    latt.construct();
    latt.print();

    Basis<Basis_t> b(&latt, nu, nd, nPho, -1);
    b.construct();
    std::cout << b << std::endl;

    double J1{1.0};
    double t1{1.0};
    double Vcu{0.0}, Vo{3.0};
    double Ucu{8.0}, Uo{4.0};
    Hamiltonian<dataType, Basis_t> H(&latt, &b, &b, true, true, 1, 1);
    // Link<dataType> J1Sq(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, J1, {{1.0,0.0,0.0}, {0.0,1.0,0.0}});
    // Link<dataType> t1Sq(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t1, {{1.0,0.0,0.0}, {0.0,1.0,0.0}});
    Link<dataType> tCuPx(LINK_TYPE::HOPPING_T, {ORBITAL::Dx2y2, ORBITAL::Px},    -t1, {{0.5, 0.0, 0.0}});
    Link<dataType> tPxCu(LINK_TYPE::HOPPING_T, {ORBITAL::Px,    ORBITAL::Dx2y2},  t1, {{0.5, 0.0, 0.0}});
    H.pushLinks({tCuPx, tPxCu}).pushV({ORBITAL::Dx2y2}, Vcu).pushV({ORBITAL::Px}, Vo).pushU({ORBITAL::Dx2y2}, Ucu).pushU({ORBITAL::Px}, Uo).transform().construct();
    H.diag();
//    H.print("Hamiltonian info");
//    std::cout << "eval: " << H.getEval() << std::endl;

    Nocc<Basis_t> nocc(&latt, &b, &b);
    nocc.construct();
    nocc.count(H.getEvec());
    std::cout << "chage distribution: " << nocc.lastCount() << std::endl;
    std::vector<std::unique_ptr<SzkOp<dataType, Basis_t>>> szk;
    for (int i = 0; i < nSite; ++i) {
        szk.push_back(std::unique_ptr<SzkOp<dataType, Basis_t>>(new SzkOp<dataType, Basis_t>(i, ORBITAL::Px, &latt, &b, &b)));
        szk[i]->construct();
    }
    std::vector<dataType> szkVal;
    for (int i = 0; i < nSite; ++i) {
        std::vector<dataType> v1(szk[i]->getnloc(), 0.0);
        szk[i]->MxV(H.getEvec(), v1.data());
        szkVal.push_back(szk[(nSite - i) % nSite]->vMv(H.getEvec(), v1.data()));
    }
    std::cout << "szk:\n" << szkVal << std::endl;
    MPI_Finalize();
    return 0;
}
