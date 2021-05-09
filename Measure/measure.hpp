#ifndef __MEASURE_H__
#define __MEASURE_H__

#include <vector>
#include <memory>
#include "Basis/Basis.hpp"
#include "Geometry/Geometry.hpp"
#include "Operator/Operators.hpp"
#include "Utils/paras.hpp"
#include "Utils/path.hpp"
#include "Solver/PARPACKSolver.hpp"
#include "Solver/Spectra.hpp"
#include "Pulse/pulse.hpp"
#include "Solver/TimeEvolver.hpp"

void init(int &workerID, int &workerNum, bool &isMaster) {
    MPI_Init(NULL, NULL);
    mpi_info(workerID, workerNum);
    isMaster = (workerID == 0);
}

/**
 * @brief sort evals and evecs in ascending order
 * 
 * @tparam T eval and evec data type
 * @param evals vector<T>
 * @param evecs  vector<T*>: pointers to evecs
 * @param groundState bool: if true return number ground state
 * @param degeneracyTol double: if abs(e1 - e2) < tol, then e1 and e2 is considered to be degenerate
 * @return int: number of states
 */
template <typename T>
int sort(std::vector<T>& evals, std::vector<T*>& evecs, bool groundState = true, double degeneracyTol = 1e-8) {
    assert_msg(evals.size() == evecs.size(), "not have same number of eval and evec");
    using eigpair = std::pair<T,T*>;
    std::vector<eigpair> pairs;
    for (size_t i = 0; i < evals.size(); ++i) {
        pairs.push_back(std::make_pair(evals[i], evecs[i]));
    }
    std::sort(pairs.begin(), pairs.end(), [](const eigpair &a, const eigpair &b) {return std::real(a.first) < std::real(b.first);});
    for (size_t i = 0; i < pairs.size(); ++i) {
        evals[i] = pairs[i].first;
        evecs[i] = pairs[i].second;
    }
    int stateNum = 0;
    if (groundState) {
        auto w0 = std::real(evals.at(0));
        for (auto &w : evals) {
            if (std::abs(std::real(w) - w0) > degeneracyTol) break;
            ++stateNum;
        }
    } else {
        stateNum = evals.size();
    }
    return stateNum;
}

enum class MEASUREMENT_TYPE {expc, corr2, spectra, time};

template <typename T>
struct Measurement {
    MEASUREMENT_TYPE type;
    OperatorBase<T>* op{nullptr};
};

// template <typename T>
// Measurement<T> keyToMeas(std::string key) {
//     tolower(key);
//     if (key == ) {

//     } else if (key)
// }

/**
 * @brief struct containing lattice, current hilbert space basis and Hamiltonian
 * 
 * @tparam T double | cdouble, Hamiltonian matrix element data type
 */
template <typename T>
struct System {
public:
    System(std::string inputDir, bool isMaster_);
    ~System( ) { }

public:
    void read(std::string inputDir);
    void construct( );
    void clear( );
    void diag( );
    bool measure(std::string key) const { return opt(measurePara, key); }

public:
    bool isMaster;
    std::unique_ptr<Geometry> latt{nullptr};
    std::unique_ptr<Basis> B{nullptr};
    std::unique_ptr<HamiltonianBase<T>> H{nullptr};

    bool isZeroTmp{true};
    int stateNum{0};
    std::vector<T> evals;
    std::vector<T*> evecs;
    int krylovDim{0};
    
    Parameters para, pulsePara, measurePara, pathPara;
    Path path;
};
template<typename T>
System<T>::System(std::string inputDir, bool isMaster_) {
    isMaster = isMaster_;
    read(inputDir);
    path = Path(&pathPara, &para, &pulsePara);
    if (isMaster) {
        path.make();
        path.print();
    }
    isZeroTmp = measure("zeroTmp");
    krylovDim = measurePara.geti("krylovDim");
    construct();
}

/**
 * @brief read input files from inputDir & construct output paths
 * 
 * @tparam T 
 * @param inputDir 
 */
template<typename T>
void System<T>::read(std::string inputDir) {
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
template<typename T>
void System<T>::construct() {
    setBasics(para, latt, B, H);
    B->construct(opt(para, "basis"), path.getBasisDir(B->getkIndex(), B->getPGIndex()));
    H->construct();
}
//TODO:System:clear
/**
 * @brief clear basis and hamiltonian matrix unused mem
 * 
 * @tparam T 
 */
template<typename T>
void System<T>::clear() {
    
}
/**
 * @brief diagonalize Hamiltonian using PARPACK. the resulting evals and evecs are sorted
 * 
 * @tparam T 
 */
template<typename T>
void System<T>::diag( ) {
    if (measure("state")) {
        evals.clear();
        evecs.clear();
        int nev = measurePara.geti("nev");
        if (H) {
            H->diag(nev);
            for (int i = 0; i < nev; ++i) {
                evals.push_back(H->getEval(i));
                evecs.push_back(H->getEvec(i));
            }
            stateNum = sort<T>(evals, evecs, isZeroTmp); 
        } else {
            std::cout<<"Warning: H is null in System.diag()!\n";
        }
    }
}

//TODO
template <typename T>
void compute(System<T> &sys, std::string key, int workerID, int workerNum, bool isSave = true) {
    assert_msg(!sys.evecs.empty(), "No states to be measured!");
    Timer timer(sys.isMaster);
    auto model = sys.B->getModel();
    for (auto &c : key) {
        c = tolower(c);
    }
    switch (model) {
        case LATTICE_MODEL::HUBBARD: {
            if (key == "conductivity") {
                Current J(sys.latt.get(), sys.B.get());
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
                auto occi = sys.B->getOcc();
                auto ki = sys.B->getkIndex();
                std::unique_ptr<Basis> Bf;
                std::unique_ptr<HamiltonianBase<T>> Hf;
                for (int kf = 0; kf < sys.latt->getSiteNum(); ++kf) {
                    for (auto spin : std::vector<SPIN> {SPIN::UP, SPIN::DOWN}) {
                        for (auto pm : std::vector<LADDER> {LADDER::PLUS, LADDER::MINUS}) {
                            auto occf = occi;
                            auto& n = (spin == SPIN::UP) ? occf.at(0) : occf.at(1);
                            auto dn = (pm == LADDER::PLUS) ? 1 : -1;
                            n += dn;
                            setbasis(sys.para, Bf, sys.latt.get(), occf.at(0), occf.at(1), kf, -1);
                            setham(sys.para, Hf, sys.latt.get(), Bf.get());
                            Bf->construct(opt(sys.para, "basis"), sys.path.getBasisDir(Bf->getkIndex(),  Bf->getPGIndex()));
                            Hf->construct();
                            for (auto& orb : sys.latt->getUnitCell()) {
                                CkOp<T> ck(sys.latt.get(), sys.B.get(), Bf.get());
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

                Nocc occ(sys.latt.get(), sys.B.get()); 
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

        case LATTICE_MODEL::HEISENBERG: {
            if (key == "sisj") {
                timer.tik();
                SSOp<T> SS(sys.latt.get(), sys.B.get());
                std::vector<std::vector<dataType>> ssvals;
                ssvals.resize(sys.stateNum);
                for (int i = 0; i < sys.latt->getSiteNum(); ++i) {
                    SS.setr(i);
                    SS.construct();
                    for (int n = 0; n < sys.stateNum; ++n) {
                        auto val = SS.vMv(sys.evecs[n], sys.evecs[n]) / dataType(sys.latt->getSiteNum());
                        // if (isMaster) std::cout<<"ss(r="<<i<<")="<<val<<std::endl;
                        ssvals[n].push_back(val);
                    }
                }
                for (int n = workerID; n < sys.stateNum; n += workerNum) {
                    save(ssvals[n].data(), (int)ssvals[n].size(), sys.path.SiSjFile + "_" + tostr(n));
                }
                timer.print("Spin Spin Correlation");

            } else if (key == "skw") {
                auto occi = sys.B->getOcc();
                for (int kf = 0; kf < sys.latt->getSiteNum(); ++kf) {
                    std::unique_ptr<Basis> Bf;
                    std::unique_ptr<HamiltonianBase<T>> Hf;
                    setbasis(sys.para, Bf, sys.latt.get(), occi.at(0), occi.at(1), kf, -1);
                    setham(sys.para, Hf, sys.latt.get(), Bf.get());
                    Bf->construct(opt(sys.para,"basis"), sys.path.getBasisDir(Bf->getkIndex(), Bf->getPGIndex()));
                    Hf->construct();
                    SzkOp<T> sk(sys.latt.get(), sys.B.get(), Bf.get());
                    sk.construct();
                    for (int n = 0; n < sys.stateNum; ++n) {
                        SPECTRASolver<dataType> spectra(Hf.get(), sys.evals[n], &sk, sys.evecs[n], sys.krylovDim);
                        spectra.compute();
                        if (sys.isMaster) spectra.save(sys.path.SkwDir + "/kf" + tostr(kf), n);
                    }
                }
            } else if (key == "raman") {
                auto J1 = sys.para.getd("J1");
                auto J2 = sys.para.getd("J2");
                for (auto channel : std::vector<std::string> {"A1", "A2", "E21", "E22"}) {
                    Hamiltonian<LATTICE_MODEL::HEISENBERG, dataType> R(sys.latt.get(), sys.B.get(), sys.B.get());
                    R.pushLinks(RamanChannel(channel, J1, J2, *(sys.latt)));
                    R.construct();
                    for (int n = 0; n < sys.stateNum; ++n) {
                        SPECTRASolver<dataType> spectra(sys.H.get(), sys.evals[n], &R, sys.evecs[n], sys.krylovDim);
                        spectra.compute();
                        if (sys.isMaster) spectra.save(sys.path.RamanDir + "/" + channel, n);
                    }
                }
            }
            break;
        }

        case LATTICE_MODEL::tJ: {
            break;
        }

        default:
            break;
    }
}

#endif // __MEASURE_H__