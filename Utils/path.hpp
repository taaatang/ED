#ifndef __PATH_H__
#define __PATH_H__

#include <iostream>
#include <string>
#include <sstream>

#include "Utils/paras.hpp"
#include "Utils/utils.hpp"

class Path {
public:
    Path(const Parameters* path, const Parameters* model, const Parameters* pulse = nullptr);
    ~Path( ) { }
    
    void make( ) const;

    void print(std::ostream& os = std::cout) const;

    std::string sep{"/"};

    std::string rootDir;
    std::string project;

    std::string lattice, particleNum, lattSymm;
    std::string symmLabel;

    std::string parameterLabel;
    std::string parameterDir;
    std::string parameterFile;

    std::string basisDir;
    std::string basisFile;
    std::string basisNormFile;

    std::string wavefuncDir;

    // Hubbard files
    std::string particleDistrDir;
    std::string sigmaDir;
    std::string AkwDir;
    std::string pumpDir;
    std::string pumpFile;

    // Heisenberg files
};

Path::Path(const Parameters* path, const Parameters* model, const Parameters* pulse) {
    rootDir = path->maps.at("rootDir");
    project = path->maps.at("project");

    lattice = model->maps.at("lattice type") + sep + model->maps.at("boundary condition");
    std::string lattSize;
    lattSize = tostr(model->mapi.at("lx")) + "x" + tostr(model->mapi.at("ly"));
    lattice += sep + lattSize;

    particleNum = tostr(model->mapi.at("nu")) + "u" + tostr(model->mapi.at("nd")) + "d";

    lattSymm = "k" + tostr(model->mapi.at("kidx"));
    auto pidx = model->mapi.at("pidx");
    if (pidx >= 0) {
        lattSymm += "p" + tostr(pidx);
    }

    symmLabel = particleNum + sep + lattSymm;

    parameterLabel = model->maps.at("parameter label");
    parameterDir = rootDir + sep + project + sep + lattice + sep + symmLabel + sep + parameterLabel;
    parameterFile = parameterDir + sep + "parameter.txt";

    basisDir = rootDir + sep + project + sep + lattice + sep + "Basis" + sep + symmLabel;
    basisFile = basisDir + sep + "basis";
    basisNormFile = basisDir + sep + "norm";

    wavefuncDir = parameterDir + sep + "wavefunc";

    switch (model->getmodel()) {
        case LATTICE_MODEL::HUBBARD:
            particleDistrDir = parameterDir + sep + "occpt";

            sigmaDir = parameterDir + sep + "sigma";

            AkwDir = parameterDir + sep + "Akw";

            if (pulse) {
                pumpDir = parameterDir + sep + "pump" + \
                sep + "pol_" + pulse->maps.at("polarization label") + \
                sep + "i_" + tostr(pulse->mapd.at("fluence")) + \
                sep + "w_" + tostr(pulse->mapd.at("frequency")) + \
                sep + "phi_" + tostr(pulse->mapd.at("phase"));

                pumpFile = pumpDir + sep + "pumpParameter.txt";
            }

            break;
        case LATTICE_MODEL::tJ:
            break;
        case LATTICE_MODEL::HEISENBERG:
            break;
        default:
            break;
    }
    
}

void Path::make( ) const { 
    if (!parameterDir.empty()) mkdir_fs(parameterDir);
    if (!basisDir.empty()) mkdir_fs(basisDir);
    if (!wavefuncDir.empty()) mkdir_fs(wavefuncDir);
    if (!particleDistrDir.empty()) mkdir_fs(particleDistrDir);
    if (!sigmaDir.empty()) mkdir_fs(sigmaDir);
    if (!AkwDir.empty()) mkdir_fs(AkwDir);
    if (!pumpDir.empty()) mkdir_fs(pumpDir);
}

void Path::print(std::ostream& os) const {
    os<<"root dir:"<<rootDir<<"\n";
    os<<"project:"<<project<<"\n";
    os<<"basis dir:"<<basisDir<<"\n";
    os<<"job data dir:"<<parameterDir<<"\n";
    if (!wavefuncDir.empty()) os<<"wavefunc dir:"<<wavefuncDir<<"\n"; 
    if (!particleDistrDir.empty()) os<<"particle distribution dir:"<<particleDistrDir<<"\n";
    if (!sigmaDir.empty()) os<<"conductivity dir:"<<sigmaDir<<"\n"; 
    if (!AkwDir.empty()) os<<"Akw dir:"<<AkwDir<<"\n"; 
    if (!pumpDir.empty()) os<<"pump dir:"<<pumpDir<<"\n"; 
}

#endif // __PATH_H__