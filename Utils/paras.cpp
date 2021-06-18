#include "Utils/paras.hpp"

Orbital stringToOrb(std::string name, int id) {
    ORBITAL orb;
    Vec3d coord;
    // tolower(name);
    if (name == "single") {
        orb = ORBITAL::SINGLE;
        coord = {0.0, 0.0, 0.0};
    } else if (name == "dx2y2") {
        orb = ORBITAL::Dx2y2;
        coord = {0.0, 0.0, 0.0};
    } else if (name == "dz2") {
        orb = ORBITAL::Dz2;
        coord = {0.0, 0.0, 0.0};
    } else if (name == "px") {
        orb = ORBITAL::Px;
        coord = {0.5, 0.0, 0.0};
    } else if (name == "py" or name == "py+") {
        orb = ORBITAL::Py;
        coord = {0.0, 0.5, 0.0};
    } else if (name == "py-") {
        orb = ORBITAL::Py;
        coord = {0.0, -0.5, 0.0};
    } else if (name == "pz" or name == "pzu") {
        orb = ORBITAL::Pzu;
        coord = {0.0, 0.0, 0.5};
    } else if (name == "pzd") {
        orb = ORBITAL::Pzd;
        coord = {0.0, 0.0, -0.5};
    } else {
        std::cout<<"Orbital: "<<name<<", is not defined!\n";
        exit(1);
    }
    return Orbital(orb, id, coord);
}

LATTICE Parameters::getlatt() const {
    auto name = get<std::string>("lattice type");
    if(name.value() == "square") {
        return LATTICE::SQUARE;
    } else if(name.value() == "triangular") {
        return LATTICE::TRIANGULAR;
    } else {
        std::cout<<"lattice type: "<<name.value()<<" is not defined!";
        exit(1);
    }
}

Parameters::Parameters(std::string inputDir, std::vector<std::string> files) {
    for(auto const& file:files){
        if(!read(inputDir+"/"+file)) {
            std::cout<<"Parameters read err.\n";
            exit(1);
        }
    }
}

LATTICE_MODEL Parameters::getmodel() const {
    auto name = get<std::string>("model name");
    if (name.value() == "Hubbard") {
        return LATTICE_MODEL::HUBBARD;
    } else if (name.value() == "Heisenberg") {
        return LATTICE_MODEL::HEISENBERG;
    } else if (name.value() == "tJ") {
        return LATTICE_MODEL::tJ;
    } else {
        std::cout << "Model name: " << name.value() << ", not defined!\n";
        exit(1);
    }
}
//TODO:Clear up
void Parameters::config(std::string configFile){
    clear();
    if(!read(configFile)) {
        std::cout << "Parameters read err.\n";
        exit(1);
    }
    assert(maps.find("rootDataPath")!=maps.end()); rootDataPath = maps["rootDataPath"];maps.erase("rootDataPath");
    assert(maps.find("project")!=maps.end()); project = maps["project"];maps.erase("project");
    assert(mapvecs.find("inputfiles")!=mapvecs.end());
    assert(maps.find("inputDir")!=maps.end()); std::string inputDir = maps["inputDir"];maps.erase("inputDir");
    assert(!mapvecs["inputfiles"].empty());
    for(auto const& file:mapvecs["inputfiles"]){
        if(!read(inputDir+"/"+file)){std::cout<<"Parameters read err.\n";exit(1);}
    }
    mapvecs.erase("inputfiles");
}


void Parameters::clear() {
    rootDataPath.clear(); 
    project.clear();
    
    mapi.clear(); 
    mapd.clear(); 
    maps.clear();

    mapvecd.clear(); 
    mapvecs.clear(); 
    maparrd.clear();
}

bool Parameters::read(const std::string& filename) {
    std::ifstream infile(filename);
    if(infile.is_open()){
        std::string line;
        while(std::getline(infile,line)){
            removeComment(line,'#');
            std::istringstream ins(line);
            std::string parsed;
            std::vector<std::string> vals;
            while(std::getline(ins,parsed,':')) {
                vals.push_back(parsed);
            }
            if(vals.size()==3) {
                std::string key = vals[0];
                std::string type = vals[1];
                std::istringstream ss(vals[2]);
                if (type == "b") {
                    bool val;
                    ss >> val;
                    mapb[key] = val;
                } else if (type == "i") {
                    int val;
                    ss >> val;
                    mapi[key] = val;
                } else if (type == "d") {
                    double val;
                    ss >> val;
                    mapd[key] = val;
                } else if (type == "s") {
                    std::string val;
                    ss >> val;
                    maps[key] = val;
                }
            } else if (vals.size() == 2) {
                std::string key = vals[0];
                std::string type = vals[1];
                std::string nextLine;
                if (type == "vecs") {
                    if (std::getline(infile,nextLine)) {
                        auto val = readVec<std::string>(nextLine);
                        mapvecs[key] = val;
                    }
                } else if (type == "vecd") {
                    if (std::getline(infile,nextLine)) {
                        auto val = readVec<double>(nextLine);
                        mapvecd[key] = val;
                    }
                } else if (type == "arrd") {
                    while (std::getline(infile,nextLine)) {
                        if(nextLine=="")break;
                        auto val = readVec<double>(nextLine);
                        maparrd[key].push_back(val);
                    }
                }
            }
        }
        infile.close();
        return true;
    } else {
        std::cout << filename << " failed to open!\n";
        return false;
    }
}

void Parameters::print(std::string filename) const {
    std::ofstream outfile(filename);
    print(outfile);
}
void Parameters::print(std::ostream& os) const {
    os << maps << mapi << mapd << mapvecs << mapvecd << maparrd;
}

void setlatt(const Parameters& para, std::unique_ptr<Geometry>& latt) {
    auto lx = para.get<int>("lx").value();
    auto ly = para.get<int>("ly").value();
    auto lattType = para.getlatt();
    bool BC = (*para.get<std::string>("boundary condition") == "periodic");
    switch (lattType) {
        case LATTICE::TRIANGULAR:
            if(ly > 0) {
                latt = std::unique_ptr<Geometry>(new TriAngLattice(lx, ly, BC));
            } else {
                latt = std::unique_ptr<Geometry>(new TriAngLattice(lx, BC));
            }
            break;
        case LATTICE::SQUARE:
            if(ly > 0) {
                latt = std::unique_ptr<Geometry>(new SquareLattice(lx, ly, BC));
            } else {
                latt = std::unique_ptr<Geometry>(new SquareLattice(lx, BC));
            }
            break;
        default:
            break;
    }

    int id = 0;
    //!what is the problem for this? tempory variable para.get<VecStr>(key) was destroyed before for-loop finished?
    // for (auto orbname : *para.get<VecStr>("orbitals")) {
    //     latt->addOrb(stringToOrb(orbname, id));
    //     id++;
    // }
    auto orbs = *para.get<VecStr>("orbitals");
    for (auto orbname : orbs) {
        latt->addOrb(stringToOrb(orbname, id));
        id++;
    }

    latt->construct();
    if (!latt->check()) {
        latt->print();
        std::cout << "Lattice does not pass check()!\n";
        exit(1);
    }
    // std::cout<<"lattice set!\n";
}

void setbasis(const Parameters& para, std::unique_ptr<Basis>& ba, Geometry* latt){
    auto kid = para.get<int>("kidx").value();
    auto pid = para.get<int>("pidx").value();
    auto nu = para.get<int>("nu").value();
    auto nd = para.get<int>("nd").value();
    ba = std::unique_ptr<Basis>(new Basis(para.getmodel(), latt, {nu, nd}, kid, pid));
    // std::cout << "basis set!\n";
}

void setbasis(std::unique_ptr<Basis>& B, LATTICE_MODEL model, Geometry* latt, int nuf, int ndf, int kf, int pf) {
    B = std::unique_ptr<Basis>(new Basis(model, latt, {nuf,ndf}, kf, pf));
    // std::cout << "basis set!\n";
}

void setham(const Parameters& para, std::unique_ptr<HamiltonianBase<dataType>>& H, Geometry* latt, Basis* B) {
    LATTICE_MODEL model = para.getmodel();
    switch (model) {
        case LATTICE_MODEL::HUBBARD:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::HUBBARD,dataType>(latt, B, B, true, 1, 1));
            break;
        case LATTICE_MODEL::tJ:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::tJ,dataType>(latt, B, B, true, 1, 1));
            break;
        case LATTICE_MODEL::HEISENBERG:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::HEISENBERG,dataType>(latt, B, B, true, 1, 1));
            break;
        default:
            std::cout << "Input Lattice Model not defined!\n";
            exit(1);
            break;
    }
    if (model == LATTICE_MODEL::HUBBARD) {
        auto unitcell = latt->getUnitCell();
        if (unitcell.size() == 1) {/* single band */
            H->pushU({ORBITAL::SINGLE}, para.get<double>("Uss").value());
            auto links = HubbardSingleBandLink(*latt);
            auto tnn = para.get<double>("tnn").value();
            auto tnnn = para.get<double>("tnnn").value();
            if (std::abs(tnn) > INFINITESIMAL) {/* nearest neighbor hopping */
                for (int id = 0; id < 2; ++id) {
                    links.at(id).setVal(links.at(id).getVal() * tnn);
                    H->pushLink(links.at(id), 0);
                }
            }
            if (std::abs(tnnn) > INFINITESIMAL) {/* next nearest neighbor hopping */
                for (int id = 2; id < 4; ++id) {
                    links.at(id).setVal(links.at(id).getVal() * tnnn);
                    H->pushLink(links.at(id), 0);
                }
            }
        } else {/* multi band */
            for (auto orb : unitcell) {
                auto ids = latt->getOrbID(orb.orb);
                auto id = ids.at(0);
                auto v = (*para.get<ArrD>("t")).at(id).at(id);
                auto u = (*para.get<ArrD>("U")).at(id).at(id);
                H->pushV({orb.orb}, v);
                H->pushU({orb.orb}, u);
            }
            auto links = HubbardMultiBandLink(*latt);
            for (auto& link : links) {
                auto orbs = link.getOrbs();
                bool validLink = true;
                VecI orbids;
                for (auto orb : orbs) {
                    auto ids = latt->getOrbID(orb);
                    if (ids.empty()) {
                        validLink = false;
                        break;
                    }
                    orbids.push_back(ids.at(0));
                }
                if (!validLink) {
                    continue;
                }
                assert(orbids.size()==2 and orbids[0]!=orbids[1]);
                std::optional<ArrD> arrdOpt = std::nullopt;
                switch (link.getLinkType()) {
                    case LINK_TYPE::HOPPING_T:
                        arrdOpt = para.get<ArrD>("t");
                        break;
                    case LINK_TYPE::HUBBARD_U:
                        arrdOpt = para.get<ArrD>("U");
                        break;
                    case LINK_TYPE::EXCHANGE_J:
                        arrdOpt = para.get<ArrD>("J");
                        break;
                    case LINK_TYPE::PAIR_HOPPING_J:
                        arrdOpt = para.get<ArrD>("J");
                        break;
                    default:
                        break;
                }
                if (arrdOpt) {
                    auto val = arrdOpt->at(orbids[0]).at(orbids[1]);
                    if (std::abs(val) > INFINITESIMAL) {
                        link.setVal(link.getVal() * val);
                        H->pushLink(link, 0);
                    }
                }
            }
        }
    } else if (model == LATTICE_MODEL::HEISENBERG) {
        auto J1 = para.get<double>("J1");
        if (std::abs(J1.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("J1", *latt);
            link.setVal(link.getVal() * J1.value_or(0.0));
            H->pushLink(link, 0);
        }
        auto J2 = para.get<double>("J2");
        if (std::abs(J2.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("J2", *latt);
            link.setVal(link.getVal() * J2.value_or(0.0));
            H->pushLink(link, 0);
        }
        auto Jk = para.get<double>("Jk");
        if (std::abs(Jk.value_or(0.0)) > INFINITESIMAL) {
            auto link = HeisenbergLink("Jk", *latt);
            link.setVal(link.getVal() * Jk.value_or(0.0));
            H->pushLink(link, 0);
        }

    }
    // std::cout<<"hamiltonian set!\n";
}

void setBasics(const Parameters& para, std::unique_ptr<Geometry>& latt, std::unique_ptr<Basis>& B, std::unique_ptr<HamiltonianBase<dataType>>& H) {
    setlatt(para, latt);
    setbasis(para, B, latt.get());
    setham(para, H, latt.get(), B.get());
}

void setpulse(const Parameters& para, Pulse& pulse) {
    auto useA = para.get<bool>("useA");
    auto freq = para.get<double>("frequency");
    auto phase = para.get<double>("phase");
    auto dt = para.get<double>("dt");
    auto duration = para.get<double>("duration");
    auto sigma = para.get<double>("sigma");
    auto fluence = para.get<double>("fluence");
    auto polv = para.get<VecD>("polarization");
    assert_msg(polv->size() == 3, "polarization should be a 3-vec!");
    Vec3d pol{polv->at(0), polv->at(1), polv->at(2)};
    pulse = Pulse(freq.value_or(0.0), sigma.value_or(50.0), dt.value_or(0.01), duration.value_or(10.0), useA.value_or(true));
    pulse.setFuncPara();
    pulse.setPol(pol);
    pulse.setPhase(phase.value_or(0.0));
    pulse.setFluence(fluence.value_or(1.0));
}

bool opt(const Parameters &para, std::string key) {
    auto cond = para.get<bool>(key, REQUIRED::NO);
    if (cond) {
        return cond.value();
    } else {
        std::cout << key << " not found in opt(para,key) function!\n";
        return false;
    }
}