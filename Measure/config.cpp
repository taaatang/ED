#include "Measure/config.hpp"

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

LATTICE getLattType(const Parameters& para) {
    auto name = para.template get<std::string>("lattice type");
    if(name.value() == "square") {
        return LATTICE::SQUARE;
    } else if(name.value() == "triangular") {
        return LATTICE::TRIANGULAR;
    } else {
        std::cout<<"lattice type: "<<name.value()<<" is not defined!";
        exit(1);
    }
}

LATTICE_MODEL getModel(const Parameters& para) {
    auto name = para.template get<std::string>("model name");
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

void setlatt(const Parameters& para, std::unique_ptr<Geometry>& latt) {
    auto lx = para.get<int>("lx").value();
    auto ly = para.get<int>("ly").value();
    auto lattType = getLattType(para);
    bool BC = *(para.template get<std::string>("boundary condition")) == "periodic";
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
    auto orbs = *(para.template get<VecStr>("orbitals"));
    for (auto orbname : orbs) {
        latt->addOrb(stringToOrb(orbname, id++));
    }

    latt->construct();
    if (!latt->check()) {
        latt->print();
        std::cout << "Lattice does not pass check()!\n";
        exit(1);
    }
}

void setpulse(const Parameters& para, Pulse& pulse) {
    auto useA = para.template get<bool>("useA");
    auto freq = para.template get<double>("frequency");
    auto phase = para.template get<double>("phase");
    auto dt = para.template get<double>("dt");
    auto duration = para.template get<double>("duration");
    auto sigma = para.template get<double>("sigma");
    auto fluence = para.template get<double>("fluence");
    auto polv = para.template get<VecD>("polarization");
    assert_msg(polv->size() == 3, "polarization should be a 3-vec!");
    Vec3d pol{polv->at(0), polv->at(1), polv->at(2)};
    pulse = Pulse(freq.value_or(0.0), sigma.value_or(50.0), dt.value_or(0.01), duration.value_or(10.0), useA.value_or(true));
    pulse.setFuncPara();
    pulse.setPol(pol);
    pulse.setPhase(phase.value_or(0.0));
    pulse.setFluence(fluence.value_or(1.0));
}