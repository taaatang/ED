//
// Created by tatang on 8/2/21.
//

#include "Basis/BasisState.hpp"
#include "Basis/Basis.hpp"
#include "Operator/localOperators.hpp"
#include "Operator/Operators.hpp"

using Basis_t = SpinBasis;

int main() {
    int workerID{0}, workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(isMaster);
    std::cout << std::setprecision(20);
    int nSite = 4;
    int nu = 2;
    int nd = 2;
    int nPho = 3;
//    SpinBasis::configure(nSite, nu, nd, nPho);
//    SpinBasis state;

//    ElectronBasis::setAllowDoubleOcc(true);
//    ElectronBasis::configure(nSite, nu, nd, nPho);
//    ElectronBasis state;

//    ElectronPhononBasis::configure(nSite, nu, nd, nPho);
//    ElectronPhononBasis state;
//    int count = 0;
//    do {
//        count++;
//        std::cout << state << std::endl;
//    } while (state.next());
//    std::cout << "dim " << state.getTotDim() << " == count " << count << std::endl;

    SquareLattice latt(nSite, 1, false);
    latt.addOrb(Orbital{ORBITAL::SINGLE});
    latt.construct();
    latt.print();
    ElectronBasis::setAllowDoubleOcc(false);
    Basis<Basis_t> b(&latt, nu, nd, nPho, -1);
    b.construct();
    std::cout << b << std::endl;
//     Test search
//        for (idx_t i = 0; i < b.getSubDim(); ++i) {
//            auto idxOpt = b.search(b.get(i));
//            assert_msg(idxOpt != std::nullopt, "not found!");
//            std::cout << "found " << *idxOpt << std::endl;
//        }

// Test localOperators
//    auto bvi = BVopt<double, ElectronPhononBasis>{b.get(760)};
//    auto bvf = Sz(2) * bvi;
//    std::cout << "initial: " << bvi->basis << ' ' << bvi->val << std::endl;
//    if (bvf) {
//        std::cout << "final  : " << bvf->basis << ' ' << bvf->val << std::endl;
//    }

    Hamiltonian<dataType, Basis_t> H(&latt, &b, &b, true, true, 1, 0);
    Link<dataType> SqLattJ1Link(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {{1.0,0.0,0.0}, {0.0,1.0,0.0}});
    H.pushLinks({SqLattJ1Link});
    H.transform();
    H.construct();
    H.diag();
    H.print("Hamiltonian info");
    std::cout << "eval: " << H.getEval() << std::endl;
    MPI_Finalize();
    return 0;
}
