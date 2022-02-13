//
// Created by Ta Tang on 2/12/22.
//
/*
 * @brief benchmark dmrg time evolution
 */
#include <iomanip>
#include "measure/measure.hpp"

using Basis_t = ElectronBasis;
int main(int argc, char *argv[]) {
    ElectronBasis::setAllowDoubleOcc(true);
    int workerID{0};
    int workerNum{1};
    bool isMaster;
    init(workerID, workerNum, isMaster);
    Timer timer(true);
    
    int L = 8;
    int nu = 4;
    int nd = 4;
    double t = -1.0;
    double U = 8.0;
    
    SquareLattice latt(L, 1, false);
    latt.addOrb({ORBITAL::SINGLE, 0, {0.0, 0.0, 0.0}}).construct();
    Basis<Basis_t> bi(&latt, nu, nd, 0);
    bi.construct();
    Basis<Basis_t> bf(&latt, nu - 1, nd, 0);
    bf.construct();

    Link<dataType> hop(LINK_TYPE::HOPPING_T, {ORBITAL::SINGLE, ORBITAL::SINGLE}, t, {{1.0, 0.0, 0.0}});

    Hamiltonian<dataType, Basis_t> Hi(&latt, &bi, &bi, true, true, 1, 1);
    Hi.pushLinks({hop}).pushU({ORBITAL::SINGLE}, U).transform().construct();
//    Hi.print("Hi Info");
    Hi.diag();
    std::cout << "ground state energy = " <<  std::setprecision(16) << Hi.getEval() << std::endl;

    Hamiltonian<dataType, Basis_t> Hf(&latt, &bf, &bf, true, true, 1, 1);
    Hf.pushLinks({hop}).pushU({ORBITAL::SINGLE}, U).transform().construct();
//    Hf.print("Hf Info");

    auto cUpi = OpR<dataType, CPlus<SPIN::UP>, Basis_t>(L/2, &latt, &bi, &bf);
    cUpi.construct();
//    std::cout << "cUp shape = " << cUp.getDim() << "  * " << cUp.getColDim() << std::endl;
//    cUp.print("cUp Info");
    int pos = 7;
    auto cUpj = OpR<dataType, CPlus<SPIN::UP>, Basis_t>(pos, &latt, &bi, &bf);
    cUpj.construct();

    int kDim = 30;
    double dt = 0.01;
    int steps = 1000;
    doubleTimeCorrelator<dataType, Basis_t>(&cUpj, &cUpi, &Hi, &Hf, Hi.getEvec(), kDim, dt, steps, "data/hubb1b/"+
                                                                                                   tostr(pos), isMaster);

    return 0;
}