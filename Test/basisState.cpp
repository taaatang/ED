//
// Created by tatang on 8/2/21.
//
#include "Basis/BasisState.hpp"

int main() {
    int nSite = 4;
    int nu = 2;
    int nd = 2;
    int nPho = 1;
//    SpinBasis::configure(nSite, nu, nd, nPho);
//    SpinBasis state;
    ElectronBasis::setAllowDoubleOcc(false);
    ElectronBasis::configure(nSite, nu, nd, nPho);
    ElectronBasis state(nSite, nu, nd);
//    ElectronPhononBasis state(nSite, nu, nd, nPho);
    int count = 0;
    do {
        count++;
        std::cout << state << std::endl;
    } while (state.next());
    std::cout << "dim " << state.getTotDim() << " == count " << count << std::endl;
    return 0;
}
