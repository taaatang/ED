#ifndef __PATH_H__
#define __PATH_H__

#include <iostream>
#include <string>
#include <sstream>

#include "Utils/paras.hpp"
#include "Utils/utils.hpp"

class Path {
public:
    Path(const Parameters* path, const Parameters* model, const Parameters* measure = nullptr, const Parameters* pulse = nullptr);
    ~Path( ) { }

    std::string sep{"/"};

    std::string rootDir;
    std::string project;

    std::string lattice, particleNum, lattSymm;
    std::string symmLabel;
    std::string parameterLabel;

    std::string parameterDir;

    std::string basisFile;
    std::string basisNormFile;

    std::string wavefuncFile;

    // Hubbard files
    std::string particleDistributionPath;
    std::string AkwFile;
    std::string pumpFile;

    // Heisenberg files
};

Path::Path(const Parameters* path, const Parameters* model, const Parameters* measure, const Parameters* pulse) {
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

    
}

#endif // __PATH_H__