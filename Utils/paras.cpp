#include "Utils/paras.hpp"

Orbital stringToOrb(std::string name, int id) {
    ORBITAL orb;
    Vec3d coord;
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
    std::string name = maps.at("lattice type");
    if(name == "square") {
        return LATTICE::SQUARE;
    } else if(name == "triangular") {
        return LATTICE::TRIANGULAR;
    } else {
        std::cout<<"lattice type: "<<name<<" is not defined!";
        exit(1);
    }
}

Parameters::Parameters(std::string inputDir, std::vector<std::string> files) {
    for(auto const& file:files){
        if(!read(inputDir+"/"+file)){std::cout<<"Parameters read err.\n";exit(1);}
    }
}

LATTICE_MODEL Parameters::getmodel() const {
    std::string name = maps.at("model name");
    if (name == "Hubbard") {
        return LATTICE_MODEL::HUBBARD;
    } else if (name == "Heisenberg") {
        return LATTICE_MODEL::HEISENBERG;
    } else if (name == "tJ") {
        return LATTICE_MODEL::tJ;
    } else {
        std::cout<<"Model name: "<<name<<", not defined!";
        exit(1);
    }
}

void Parameters::config(std::string configFile){
    clear();
    if(!read(configFile)){std::cout<<"Parameters read err.\n";exit(1);}
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


void Parameters::clear(){
    rootDataPath.clear(); project.clear();
    mapi.clear(); mapd.clear(); maps.clear();
    mapvecd.clear(); mapvecs.clear(); maparrd.clear();
}

bool Parameters::read(const std::string& filename){
    std::ifstream infile(filename);
    if(infile.is_open()){
        std::string line;
        while(std::getline(infile,line)){
            removeComment(line,'#');
            std::istringstream ins(line);
            std::string parsed;
            std::vector<std::string> vals;
            while(std::getline(ins,parsed,':')){
                vals.push_back(parsed);
            }
            if(vals.size()==3){
                std::string key = vals[0];
                std::string type = vals[1];
                std::istringstream ss(vals[2]);
                if(type=="i"){int val;ss>>val;mapi[key]=val;}
                else if(type=="d"){double val;ss>>val;mapd[key]=val;}
                else if(type=="s"){std::string val;ss>>val;maps[key]=val;}
            }else if(vals.size()==2){
                std::string key = vals[0];
                std::string type = vals[1];
                std::string nextLine;
                if(type=="vecs"){
                    if(std::getline(infile,nextLine)){
                        auto val = readVec<std::string>(nextLine);
                        mapvecs[key] = val;
                    }
                }else if(type=="vecd"){
                    if(std::getline(infile,nextLine)){
                        auto val = readVec<double>(nextLine);
                        mapvecd[key] = val;
                    }
                }else if(type=="arrd"){
                    while(std::getline(infile,nextLine)){
                        if(nextLine=="")break;
                        auto val = readVec<double>(nextLine);
                        maparrd[key].push_back(val);
                    }
                }
            }
        }
        infile.close();
        return true;
    }else{
        std::cout<<filename<<" failed to open!\n";
        return false;
    }
}

void Parameters::print(std::string filename) const {
    std::ofstream outfile(filename);
    print(outfile);
}
void Parameters::print(std::ostream& os) const {
    os<<maps<<mapi<<mapd<<mapvecs<<mapvecd<<maparrd;
}

void setlatt(const Parameters& para, std::unique_ptr<Geometry>& latt){
    auto lx = para.mapi.at("lx");
    auto ly = para.mapi.at("ly");
    auto lattType = para.getlatt();
    bool BC = (para.maps.at("boundary condition") == "periodic");
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
    for (auto orbname : para.mapvecs.at("orbitals")) {
        latt->addOrb(stringToOrb(orbname, id));
        id++;
    }

    latt->construct();
    if (!latt->check()) {
        latt->print();
        std::cout<<"Lattice does not pass check()!\n";
        exit(1);
    }
    // std::cout<<"lattice set!\n";
}

void setbasis(const Parameters& para, std::unique_ptr<Basis>& ba, Geometry* latt){
    int kid = para.mapi.at("kidx");
    int pid = para.mapi.at("pidx");
    int nu = para.mapi.at("nu");
    int nd = para.mapi.at("nd");
    ba = std::unique_ptr<Basis>(new Basis(para.getmodel(), latt, {nu,nd}, kid, pid));
    // std::cout<<"basis set!\n";
}

void setbasis(const Parameters& para, std::unique_ptr<Basis>& ba, Geometry* latt, int nuf, int ndf, int kf, int pf) {
    ba = std::unique_ptr<Basis>(new Basis(para.getmodel(), latt, {nuf,ndf}, kf, pf));
    // std::cout<<"basis set!\n";
}

void setham(const Parameters& para, std::unique_ptr<HamiltonianBase<dataType>>& H, Geometry* latt, Basis* B) {
    LATTICE_MODEL model = para.getmodel();
    switch (model) {
        case LATTICE_MODEL::HUBBARD:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::HUBBARD,dataType>(latt, B, B, 1, 1));
            break;
        case LATTICE_MODEL::tJ:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::tJ,dataType>(latt, B, B, 1, 1));
            break;
        case LATTICE_MODEL::HEISENBERG:
            H = std::unique_ptr<HamiltonianBase<dataType>>(new Hamiltonian<LATTICE_MODEL::HEISENBERG,dataType>(latt, B, B, 1, 1));
            break;
        default:
            std::cout<<"Input Lattice Model not defined!\n";
            exit(1);
            break;
    }
    if (model == LATTICE_MODEL::HUBBARD) {
        auto unitcell = latt->getUnitCell();
        for (auto orb : unitcell) {
            auto ids = latt->getOrbID(orb.orb);
            auto id = ids.at(0);
            auto v = para.maparrd.at("t").at(id).at(id);
            auto u = para.maparrd.at("U").at(id).at(id);
            H->pushV({orb.orb}, v);
            H->pushU({orb.orb}, u);
        }
        auto links = HubbardLink();
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
            const std::vector<std::vector<double>>* arrdPtr = nullptr;
            switch (link.getLinkType()) {
                case LINK_TYPE::HOPPING_T:
                    arrdPtr = &para.maparrd.at("t");
                    break;
                case LINK_TYPE::HUBBARD_U:
                    arrdPtr = &para.maparrd.at("U");
                    break;
                case LINK_TYPE::EXCHANGE_J:
                    arrdPtr = &para.maparrd.at("J");
                    break;
                case LINK_TYPE::PAIR_HOPPING_J:
                    arrdPtr = &para.maparrd.at("J");
                    break;
                default:
                    break;
            }
            if (arrdPtr) {
                auto val = arrdPtr->at(orbids[0]).at(orbids[1]);
                if (std::abs(val) > INFINITESIMAL) {
                    link.setVal(link.getVal() * val);
                    H->pushLink(link, 0);
                }
            }
        }
    } else if (model == LATTICE_MODEL::HEISENBERG) {
        auto links = HeisenbergLink();

    }

    // std::cout<<"hamiltonian set!\n";
}

void setBasics(const Parameters& para, std::unique_ptr<Geometry>& latt, std::unique_ptr<Basis>& B, std::unique_ptr<HamiltonianBase<dataType>>& H, int workerID) {
    setlatt(para, latt);
    setbasis(para, B, latt.get());
    setham(para, H, latt.get(), B.get());

    Timer timer;
    timer.tik();
    B->gen();
    timer.tok();
    if (workerID == MPI_MASTER) {
        timer.print("Basis construction");
        B->print();
    }

    timer.tik();
    H->construct();
    timer.tok();
    if (workerID == MPI_MASTER) {
        timer.print("Hamiltonian construction");
        H->print("Hamiltonian from worker " + tostr(workerID));
    }
    
}

void setpulse(const Parameters& para, Pulse& pulse) {
    auto freq = para.mapd.at("frequency");
    auto phase = para.mapd.at("phase");
    auto dt = para.mapd.at("dt");
    auto numSteps = para.mapi.at("numSteps");
    auto width = para.mapd.at("width");
    auto fluence = para.mapd.at("fluence");
    auto polv = para.mapvecd.at("polarization");
    assert_msg(polv.size() == 3, "polarization should be a 3-vec!");
    Vec3d pol{polv[0], polv[1], polv[2]};
    pulse = Pulse(freq, width, dt, numSteps);
    pulse.setFuncPara();
    pulse.setPol(pol);
    pulse.setFluence(fluence);
    pulse.setPhase(phase);
}