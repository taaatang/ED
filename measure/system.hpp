#pragma once

#include "measure/measure.hpp"
#include "measure/config.hpp"

/**
 * @brief contain lattice, basis and Hamiltonian; parameters and data saving path
 * 
 * @tparam T operator element type
 * @tparam B basis type
 */
template <typename T, IsBasisType B>
class System {
public:
    System(std::string inputDir, bool isMaster_);

public:
    void construct( );

    void clear( );

    void diag( );

    bool measure(std::string key) const { return opt(measurePara, key); }

    template <typename U>
    U getMpara(std::string key) const { return measurePara.template get<U>(key).value(); }

    template <typename U>
    U getPara(std::string key) const { return para.template get<U>(key).value(); }

private:
    void read(std::string inputDir);

public:
    bool isMaster;

    Timer timer;

    std::unique_ptr<Geometry> latt{nullptr};

    std::unique_ptr<Basis<B>> basis{nullptr};

    std::unique_ptr<Hamiltonian<T, B>> H{nullptr};

    std::unique_ptr<Basis<B>> bf{nullptr};

    std::unique_ptr<Hamiltonian<T, B>> Hf{nullptr};

    bool isZeroTmp{true};

    int stateNum{0};

    std::vector<T> evals;

    std::vector<T*> evecs;

    int krylovDim{0};
    
    Parameters para;

    Parameters pulsePara;

    Parameters measurePara;

    Parameters pathPara;

    Path path;

private:
    bool isDiag{false};

    bool isConstructed{false};
};

template<typename T, IsBasisType B>
System<T, B>::System(std::string inputDir, bool isMaster_): isMaster(isMaster_), timer(isMaster_) {
    read(inputDir);
    path = Path(&pathPara, &para, &pulsePara);
    if (isMaster) {
        path.make(measurePara);
        path.print();
        para.print(path.parameterFile);
    }
    isZeroTmp = measure("zeroTmp");
    auto krylovOpt = measurePara.template get<int>("krylovDim", REQUIRED::NO);
    krylovDim = krylovOpt.value_or(400);
    setlatt(para, latt);
    //FIXME: Sparse Matrix setDm init diagValist with B->locDim. Error: 1.B not constructed, locDim is totDim; 2. B constructed, locDim is subDim, not H.locDim.
    // setBasics(para, latt, B, H);
}

/**
 * @brief read input files from inputDir & construct output paths
 * 
 * @tparam T 
 * @param inputDir 
 */
template<typename T, IsBasisType B>
void System<T, B>::read(std::string inputDir) {
    pathPara = Parameters(inputDir, {"path.txt"});
    para = Parameters(inputDir, {"lattice.txt", "hamiltonian.txt"});
    pulsePara = Parameters(inputDir, {"pulse.txt"});
    measurePara = Parameters(inputDir, {"measure.txt"});
}

/**
 * @brief construct basis and Hamiltonian matrix
 * 
 * @tparam T 
 */
template<typename T, IsBasisType B>
void System<T, B>::construct() {
    if (isConstructed) {
        return;
    }
    timer.tik();
    setbasis(para, basis, latt.get());
    basis->construct(opt(para, "basis"), path.getBasisDir(B->getKIdx, B->getPIdx()));
    timer.print("Basis construction");
    timer.tik();
    setham(para, H, latt.get(), B.get());
    H->construct();
    timer.print("Hamiltonian construction");
    if (isMaster) {
        B->print();
        H->print("Hamiltonian from master worker");
    }
    isConstructed = true;
}
//TODO:System:clear
/**
 * @brief clear basis and hamiltonian matrix unused mem
 * 
 * @tparam T 
 */
template<typename T, IsBasisType B>
void System<T, B>::clear() {
    
}
/**
 * @brief diagonalize Hamiltonian using PARPACK. the resulting evals and evecs are sorted
 * 
 * @tparam T 
 */
template<typename T, IsBasisType B>
void System<T, B>::diag( ) {
    if (isDiag) {
        return;
    }
    if (measure("state")) {
        evals.clear();
        evecs.clear();
        auto nevopt = measurePara.template get<int>("nev", REQUIRED::NO);
        int nev = nevopt.value_or(1);
        if (H) {
            H->diag(nev);
            isDiag = true;
            for (int i = 0; i < nev; ++i) {
                evals.push_back(H->getEval(i));
                evecs.push_back(H->getEvec(i));
            }
            stateNum = sort<T>(evals, evecs, isZeroTmp); 
            if (isMaster) {
                std::cout << "sorted evals: " << evals << '\n';
                save(evals.data(), int(evals.size()), path.evalFile);
            }
            //? clear sparse matrix content
            // H->clear();
        } else {
            std::cout<<"Warning: H is null in System.diag()!" << std::endl;
        }
    }
}

//TODO
template <typename T, IsBasisType B>
void compute(System<T, B> &sys, std::string key, int workerID, int workerNum, bool isSave = true) {
    Timer timer(sys.isMaster);
    for (auto &c : key) {
        c = tolower(c);
    }
    if (key == "basis") {
        auto kmin_opt = sys.measurePara.template get<int>("kmin");
        auto kmax_opt = sys.measurePara.template get<int>("kmax");
        auto pmin_opt = sys.measurePara.template get<int>("pmin");
        auto pmax_opt = sys.measurePara.template get<int>("pmax");
        auto kmin = kmin_opt.value_or(0);
        auto kmax = kmax_opt.value_or(0);
        auto pmin = pmin_opt.value_or(0);
        auto pmax = pmax_opt.value_or(0);
        int count{0};
        for (auto k = kmin; k <= kmax; ++k) {
            for (auto p = pmin; p <= pmax; ++p) {
                if (count % workerNum == workerID) {
                    Basis b(sys.para.getmodel(), sys.latt.get(), {sys.para.template get<int>("nu").value(), sys.para.template get<int>("nd").value()}, k, p);
                    b.construct();
                    b.save(sys.path.getBasisDir(k, p));
                }
                count++;
            }
        }
        return;
    } 

    sys.construct();
    sys.diag();
    assert_msg(!sys.evecs.empty(), "No states to be measured!");
    auto model = sys.B->getModel();
    switch (model) { 
        case MODEL::HUBBARD: {
            assert_msg(sys.B->getPGIndex() == -1, "HUBBARD model calculation only defined with translation symm!");
            if (key == "conductivity") {
                if (sys.H->isEmpty()) {
                    sys.H->construct();
                }
                Current J(sys.latt.get(), sys.B.get(), sys.B.get(), true, false);
                J.pushLinks(HubbardLink(*(sys.latt)));
                auto ds = std::vector<std::string> {"x", "y", "z"};
                for (auto d : ds) {
                    J.setDirection(d);
                    // J.print();
                    J.construct();
                    for (int n = 0; n < sys.stateNum; ++n){
                        SPECTRASolver<cdouble> spectra(sys.H.get(), sys.evals.at(n), &J, sys.evecs.at(n), sys.krylovDim);
                        spectra.compute();
                        if (sys.isMaster && isSave) {
                            spectra.save(sys.path.sigmaDir + "/" + d, n);
                        }
                    }
                }
            } else if (key == "akw") {
                //TODO: fix Cp/m with PG symm
                auto occi = sys.B->getOcc();
                auto ki = sys.B->getkIndex();
                std::unique_ptr<Basis> Bf;
                std::unique_ptr<HamiltonianBase<T>> Hf;
                for (int kf = 0; kf < sys.latt->getUnitCellNum(); ++kf) {
                    for (auto spin : std::vector<SPIN> {SPIN::UP, SPIN::DOWN}) {
                        for (auto pm : std::vector<LADDER> {LADDER::PLUS, LADDER::MINUS}) {
                            auto occf = occi;
                            auto& n = (spin == SPIN::UP) ? occf.at(0) : occf.at(1);
                            auto dn = (pm == LADDER::PLUS) ? 1 : -1;
                            n += dn;
                            setbasis(Bf, sys.para.getmodel(), sys.latt.get(), occf.at(0), occf.at(1), kf, -1);
                            setham(sys.para, Hf, sys.latt.get(), Bf.get());
                            Bf->construct(opt(sys.para, "basis"), sys.path.getBasisDir(Bf->getkIndex(),  Bf->getPGIndex()));
                            Hf->construct();
                            for (auto& orb : sys.latt->getUnitCell()) {
                                CkOp<T> ck(sys.latt.get(), sys.B.get(), Bf.get(), true, false);
                                ck.set(pm, spin, orb);
                                ck.construct();
                                for (int n = 0; n < sys.stateNum; ++n) {
                                    SPECTRASolver<T> spectra(Hf.get(), sys.evals.at(n), &ck, sys.evecs.at(n), sys.krylovDim);
                                    spectra.compute();
                                    if (sys.isMaster && isSave) {
                                        std::ostringstream os;
                                        os<<"/ki"<<tostr(ki)<<"/kf"<<tostr(kf)<<"/orb"<<orb.orbid<<"/"<<spin<<"_"<<pm;
                                        spectra.save(sys.path.AkwDir + os.str(), n);
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (key == "pump") {
                Pulse pulse;
                setpulse(sys.pulsePara, pulse);
                if (sys.isMaster) pulse.print();
                sys.H->setPeierls(&pulse);
                // H->printPeierls();
                sys.H->construct();
                // if (isMaster) {
                //     timer.print("Set H(t)");
                //     H->print("H(t) from worker " + tostr(workerID));
                // }

                int krylovDim = 15;
                TimeEvolver<cdouble> Tevol(sys.evecs.at(0), sys.H.get(), krylovDim);

                Nocc occ(sys.latt.get(), sys.B.get(), sys.B.get()); 
                occ.construct();

                timer.tik();
                if (sys.isMaster) {
                    pulse.progressBar(100);
                }
                while (sys.H->next()) {
                    Tevol.evolve(pulse.getdt());
                    occ.count(Tevol.getVec());
                    if (sys.isMaster) {
                        pulse.progress();
                    }
                }
                timer.tok();
                if (sys.isMaster) {
                    timer.print("Timer evolution");
                }
                if (sys.isMaster && isSave) {
                    occ.save(sys.path.pumpDir);
                }
            }
            break;
        }

        case MODEL::HEISENBERG: {
            if (key == "sisj") {
                // SSOp<T> SS(sys.latt.get(), sys.B.get(), sys.B.get(), false, false);
                std::vector<std::vector<dataType>> ssvals;
                ssvals.resize(sys.stateNum);
                for (int i = 0; i < sys.latt->getUnitCellNum(); ++i) {
                    timer.tik();
                    Hamiltonian<MODEL::HEISENBERG, cdouble> SS(sys.latt.get(), sys.B.get(), sys.B.get(), true, false);
                    Link<cdouble> J(LINK_TYPE::SUPER_EXCHANGE_J, {ORBITAL::SINGLE, ORBITAL::SINGLE}, 1.0, {sys.latt->getSiteR(i)});
                    SS.pushLink(J, 0);
                    // SS.setr(i);
                    SS.transform();
                    SS.construct();
                    timer.print("S(r=" + tostr(i) + ") construction");
                    for (int n = 0; n < sys.stateNum; ++n) {
                        timer.tik();
                        auto val = SS.vMv(sys.evecs[n], sys.evecs[n]) / dataType(sys.latt->getUnitCellNum());
                        timer.print("S(r=" + tostr(i) + ", n=" + tostr(n) + ") = " + tostr(std::real(val)));
                        ssvals[n].push_back(val);
                    }
                }
                for (int n = workerID; n < sys.stateNum; n += workerNum) {
                    save(ssvals[n].data(), (int)ssvals[n].size(), sys.path.SiSjFile + "_" + tostr(n));
                }
            } else if (key == "skw") {
                //TODO: fix PG symm for skw
                // assert_msg(sys.B->getPGIndex() == -1, "skw only defined for traslation symm!");
                auto occi = sys.B->getOcc();
                for (int kf = 0; kf < sys.latt->getUnitCellNum(); ++kf) {
                    if (sys.isMaster) printLine(40, '-');
                    std::unique_ptr<Basis> Bf;
                    std::unique_ptr<HamiltonianBase<T>> Hf;
                    timer.tik();
                    setbasis(Bf, sys.para.getmodel(), sys.latt.get(), occi.at(0), occi.at(1), kf, -1);
                    Bf->construct(opt(sys.para,"basis"), sys.path.getBasisDir(Bf->getkIndex(), Bf->getPGIndex()));
                    timer.print("Bf(k=" + tostr(kf) + ") construction");
                    timer.tik();
                    setham(sys.para, Hf, sys.latt.get(), Bf.get());    
                    Hf->construct();
                    timer.print("Hf construction");
                    if (sys.isMaster) {
                        Bf->print();
                        Hf->print("Hf from master worker");
                    }
                    timer.tik();
                    SzkOp<T> sk(sys.latt.get(), sys.B.get(), Bf.get(), false, false);
                    sk.construct();
                    // release Bf resources
                    Bf.reset(nullptr);
                    timer.print("Sq(kf=" + tostr(kf) + ") construction");
                    if (sys.measure("lanczos")) {
                        if (sys.isMaster) printLine(20, '-');
                        timer.tik();
                        for (int n = 0; n < sys.stateNum; ++n) {
                            SPECTRASolver<dataType> spectra(Hf.get(), sys.evals[n], &sk, sys.evecs[n], sys.krylovDim);
                            spectra.compute();
                            if (sys.isMaster) spectra.save(sys.path.SkwDir + "/Lanczos/kf" + tostr(kf), n);
                        }
                        timer.print("Sqw Lanczos");
                    }
                    if (sys.measure("bicg")) {
                        if (sys.isMaster) printLine(20, '-');
                        timer.tik();
                        for (int n = 0; n < sys.stateNum; ++n) {
                            SPECTRASolverBiCGSTAB spectra(Hf.get(), &sk, sys.evecs[n], std::real(sys.evals[n]), n, sys.template getMpara<int>("iterMax"), sys.template getMpara<double>("wmin"), sys.template getMpara<double>("wmax"), sys.template getMpara<double>("dw"), sys.template getMpara<double>("epsilon"));
                            spectra.compute(sys.path.SkwDir + "/BiCGSTAB/kf" + tostr(kf));
                        }
                        timer.print("Sqw BiCGSTAB");
                    }
                }
            } else if (key == "raman") {
                auto J1 = sys.template getPara<double>("J1");
                auto J2 = sys.template getPara<double>("J2");
                auto ki = sys.template getPara<int>("kidx");
                auto occi = sys.B->getOcc();

                auto channels = sys.template getMpara<VecStr>("RamanChannels");
                auto pfmin = sys.template getMpara<int>("pfmin");
                auto pfmax = sys.template getMpara<int>("pfmax");
                
                for (int p = pfmin; p <= pfmax; ++p) {
                    if (sys.isMaster) printLine(40, '-');
                    std::unique_ptr<Basis> Bf;
                    std::unique_ptr<HamiltonianBase<T>> Hf;
                    timer.tik();
                    setbasis(Bf, sys.para.getmodel(), sys.latt.get(), occi.at(0), occi.at(1), ki, p);
                    Bf->construct(opt(sys.para,"basis"), sys.path.getBasisDir(Bf->getkIndex(), Bf->getPGIndex()));
                    timer.print("Bf(p=" + tostr(p) + ") construction");
                    timer.tik();
                    setham(sys.para, Hf, sys.latt.get(), Bf.get());
                    Hf->transform();
                    Hf->construct();
                    timer.print("Hf(p=" + tostr(p) + ") construction");
                    if (sys.isMaster) {
                        Bf->print();
                        Hf->print("Hf from master worker");
                    }
                    // bool commuteWithSymm = (channel == "A1" || sys.B->getPGIndex() == -1);
                    for (auto channel : channels) {
                        if (sys.isMaster) printLine(20, '-');
                        timer.tik();                       
                        Hamiltonian<MODEL::HEISENBERG, dataType> R(sys.latt.get(), sys.B.get(), Bf.get(), true, false);
                        R.pushLinks(RamanChannel(channel, J1, J2, *(sys.latt)));
                        R.transform();
                        R.construct();
                        timer.print("Raman Op channel " + channel + " construction");
                        auto label = channel;
                        if (p >= 0) {
                            label = channel + "_p" + tostr(p);
                        }
                        if (sys.measure("lanczos")) {
                            if (sys.isMaster) printLine(10, '-');
                            timer.tik();
                            for (int n = 0; n < sys.stateNum; ++n) {
                                SPECTRASolver<dataType> spectra(Hf.get(), sys.evals[n], &R, sys.evecs[n], sys.krylovDim);
                                spectra.compute();
                                if (sys.isMaster) spectra.save(sys.path.RamanDir + "/Lanczos/" + label, n);
                            }
                            timer.print("Raman Lanczos");
                        } 
                        if (sys.measure("bicg")) {
                            if (sys.isMaster) printLine(10, '-');
                            timer.tik();
                            for (int n = 0; n < sys.stateNum; ++n) {
                                SPECTRASolverBiCGSTAB spectra(Hf.get(), &R, sys.evecs[n], std::real(sys.evals[n]), n, sys.template getMpara<int>("iterMax"), sys.template getMpara<double>("wmin"), sys.template getMpara<double>("wmax"), sys.template getMpara<double>("dw"), sys.template getMpara<double>("epsilon"));
                                spectra.compute(sys.path.RamanDir + "/BiCGSTAB/" + label);
                            }
                            timer.print("Raman BiCGSTAB");
                        }
                    }
                }
            } else if (key == "ramanfl") {
                //TODO: fix for PG symm
                auto J1 = sys.para.template get<double>("J1");
                auto J2 = sys.para.template get<double>("J2");
                Vec3d plzX{1.0, 0.0, 0.0}, plzY{-1.0, 2.0, 0.0};
                std::vector<Vec3d> plz{plzX, plzY};
                std::vector<std::string> plzLabel{"x", "y"};
                for(int i = 0; i < 1; ++i) {
                    for(int j = 0; j < 2; ++j) {
                        RamanOp<dataType> R(sys.latt.get(), sys.B.get(), sys.B.get(), true, false);
                        if (std::abs(J1.value_or(0.0)) > INFINITESIMAL) {
                            auto link = triangleHeisenbergLink("J1", *sys.latt.get());
                            link.setVal(link.getVal() * J1.value_or(0.0));
                            R.pushLink(link, 0);
                        }
                        if (std::abs(J2.value_or(0.0)) > INFINITESIMAL) {
                            auto link = triangleHeisenbergLink("J2", *sys.latt.get());
                            link.setVal(link.getVal() * J2.value_or(0.0));
                            R.pushLink(link, 0);
                        }
                        R.setplz(plz[i],plz[j]);
                        R.construct();
                        for (int n = 0; n < sys.stateNum; ++n) {
                            SPECTRASolver<dataType> spectra(sys.H.get(), sys.evals[n], &R, sys.evecs[n], sys.krylovDim);
                            spectra.compute();
                            if (sys.isMaster) spectra.save(sys.path.RamanDir + "/" + plzLabel[i] + plzLabel[j], n);
                        }
                    }
                }
            }
            break;
        }

        case MODEL::tJ: {
            break;
        }

        default:
            break;
    }
}