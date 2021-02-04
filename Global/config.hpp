#include "../Global/global.hpp"
#include "../Utils/paras.hpp"
#include "../Geometry/Geometry.hpp"
#include "../Operator/OperatorBasics.hpp"
#include "../Operator/links.hpp"
#include <memory>

void setpath(const Parameters&);
void setlatt(const Parameters&, std::unique_ptr<Geometry>& latt);
void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*);
void setbasis(const Parameters&, std::unique_ptr<Basis>&, Geometry*, int nuf, int ndf, int kf, int pf);
void setham(const Parameters&, OperatorBase<dataType>& H);
void setmeasure(const Parameters&);

std::string configFile = "../Input/config.txt";
Parameters para(configFile);
std::unique_ptr<Geometry> latt;
std::unique_ptr<Basis> Bi;
std::unique_ptr<Basis> Bf;


void setpath(Parameters& para){
    
}

void setlatt(Parameters& para, std::unique_ptr<Geometry>& latt){
    auto lx = para.mapi.at("lx");
    auto ly = para.mapi.at("ly");
    auto lattType = getlatt();
    bool BC = (maps.at("boundary condition") == "periodic");
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
    for (auto orbname : mapvecs.at("orbitals")) {
        latt->addOrb(stringToOrb(orbname, id));
        id++;
    }

    latt->construct();
}

void setbasis(Parameters& para, std::unique_ptr<Basis>& ba, Geometry* latt){
    int kid = mapi["kid"];
    int pid = mapi["pid"];
    int nu = mapi["nu"];
    int nd = mapi["nd"];
    ba = std::unique_ptr<Basis>(new Basis(para.getmodel(), latt, {nu,nd}, kid, pid));
}

void setbasis(Parameters& para, std::unique_ptr<Basis>& ba, Geometry* latt, int nuf, int ndf, int kf, int pf) {
    ba = std::unique_ptr<Basis>(new Basis(para.getmodel(), latt, {nuf,ndf}, kf, pf));
}

void setham(Parameters& para, OperatorBase<dataType>& H, Geometry* latt, Basis* B) {
    LATTICE_MODEL model = para.getmodel();
    switch (model) {
        case LATTICE_MODEL::HUBBARD:
            H = std::unique_ptr<OperatorBase>(new Hamiltonian<HUBBARD,dataType>(latt, B, B, 1, 1));
            break;
        case LATTICE_MODEL::tJ:
            H = std::unique_ptr<OperatorBase>(new Hamiltonian<tJ,dataType>(latt, B, B, 1, 1));
            break;
        case LATTICE_MODEL::HEISENBERG:
            H = std::unique_ptr<OperatorBase>(new Hamiltonian<HEISENBERG,dataType>(latt, B, B, 1, 1));
            break;
        default:
            std::cout<<"Input Lattice Model not defined!\n";
            exit(1);
            break;
    }
    if (model == LATTICE_MODEL::HUBBARD) {
        auto unitcell = latt->getUnitCell();
        for (auto orb : unitcell) {
            auto ids = latt->getOrbID(orb);
            auto id = ids.at(0);
            auto v = para.maparrd.at("t").at(id).at(id);
            auto u = para.maparrd.at("U").at(id).at(id);
            H->pushV({orb}, v);
            H->pushU({orb}, u);
        }
        auto links = HubbardLink();
        for (auto& link : links) {
            auto orbs = link.getOrbs();
            bool validLink = true;
            VecI orbids;
            for (auto orb : orbs) {
                auto ids = latt->getOrbID();
                if (ids.empty()) {
                    validLink = false;
                    break;
                }
                orbids.push_back(ids.at(0));
            }
            if (!validLink) {
                break;
            }
            assert(orbids.size()==2 and orbids[0]!=orbids[1]);

            std::vector<std::vector<double>>* arrdPtr = nullptr;
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
                link.setVal(link.getVal() * val);
                H.pushLink(link, 0);
            }
        }
    } else if (model == LATTICE_MODEL::HEISENBERG) {
        auto links = HeisenbergLink();

    }
}