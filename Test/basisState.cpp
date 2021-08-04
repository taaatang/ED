//
// Created by tatang on 8/2/21.
//

#include "Basis/BasisState.hpp"
#include "Basis/Basis.hpp"
#include "Operator/localOperators.hpp"

int main() {
    int nSite = 4;
    int nu = 1;
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

    SquareLattice latt(nSite, 1);
    latt.addOrb(ORBITAL::SINGLE);
    latt.construct();
    latt.print();
    ElectronBasis::setAllowDoubleOcc(false);
    Basis<ElectronPhononBasis> b(&latt, nu, nd, nPho, 3);
    b.construct();
    std::cout << b << std::endl;

//    CPlus<SPIN::UP> cpu;
//    CPlus<SPIN::DOWN> cpd;
//    CMinus<SPIN::UP> cmu;
//    CMinus<SPIN::DOWN> cmd;
//    NCharge<SPIN::UP> nelu;
//    NCharge<SPIN::DOWN> neld;

    auto bvi = BVopt<ElectronPhononBasis>{b.get(760)};
    auto bvf = Sz(2) * bvi;
    std::cout << "initial: " << bvi->basis << ' ' << bvi->val << std::endl;
    if (bvf) {
        std::cout << "final  : " << bvf->basis << ' ' << bvf->val << std::endl;
    }
    // test search
//    for (idx_t i = 0; i < b.getSubDim(); ++i) {
//        auto idxOpt = b.search(b.get(i));
//        assert_msg(idxOpt != std::nullopt, "not found!");
//        std::cout << "found " << *idxOpt << std::endl;
//    }


    return 0;
}
