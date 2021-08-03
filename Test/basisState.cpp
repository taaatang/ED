//
// Created by tatang on 8/2/21.
//

#include "Basis/BasisState.hpp"
#include "Basis/BasisBase.hpp"

int main() {
    int nSite = 4;
    int nu = 2;
    int nd = 2;
    int nPho = 1;
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

    // test search
    for (idx_t i = 0; i < b.getSubDim(); ++i) {
        auto idxOpt = b.search(b.get(i));
        assert_msg(idxOpt != std::nullopt, "not found!");
        std::cout << "found " << *idxOpt << std::endl;
    }
    return 0;
}
