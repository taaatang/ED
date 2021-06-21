#ifndef __PATH_H__
#define __PATH_H__

#include <iostream>
#include <string>
#include <sstream>

#include "Utils/paras.hpp"
#include "Utils/io.hpp"

class Path {
public:
    Path( ) { }
    Path(const Parameters* path, const Parameters* model, const Parameters* pulse = nullptr);
    ~Path( ) { }
    
    void make(const Parameters& mPara) const;

    void print(std::ostream& os = std::cout, bool brief = true) const;

    std::string getBasisDir(int kidx, int pidx = -1) const;

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
    std::string evalFile;

    // Hubbard files
    std::string particleDistrDir;
    std::string sigmaDir;
    std::string AkwDir;
    std::string pumpDir;
    std::string pumpFile;

    // Heisenberg files
    std::string SiSjDir;
    std::string SiSjFile;
    std::string SkwDir;
    std::string RamanDir;
};

Path::Path(const Parameters* path, const Parameters* model, const Parameters* pulse) {
    rootDir = path->maps.at("rootDir");
    project = path->maps.at("project");

    lattice = model->maps.at("lattice type") + sep + model->maps.at("boundary condition");
    std::string lattSize;
    int lx = model->mapi.at("lx");
    int ly = model->mapi.at("ly");
    lattSize = tostr(lx);
    if (ly > 0) lattSize += "x" + tostr(ly);
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
    evalFile = wavefuncDir + sep + "eval";

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
            SiSjDir = parameterDir + sep + "sisj";
            SiSjFile = SiSjDir + sep + "sisj";
            SkwDir = parameterDir + sep + "skw";
            RamanDir = parameterDir + sep + "raman";
            break;
        default:
            break;
    }
    
}

std::string Path::getBasisDir(int kidx, int pidx) const {
    std::string lattSymm_ = "k" + tostr(kidx);
    if (pidx >= 0) {
        lattSymm_ += "p" + tostr(pidx);
    }

    std::string symmLabel_ = particleNum + sep + lattSymm_; 
    std::string bDir = rootDir + sep + project + sep + lattice + sep + "Basis" + sep + symmLabel_;
    return bDir;
}

void Path::make(const Parameters& mPara) const { 
    if (!parameterDir.empty()) mkdir_fs(parameterDir);
    if (!basisDir.empty() && opt(mPara, "basis")) mkdir_fs(basisDir);
    if (!wavefuncDir.empty() && opt(mPara, "state")) mkdir_fs(wavefuncDir);
    // hubbard dirs
    if (!particleDistrDir.empty()) mkdir_fs(particleDistrDir);
    if (!sigmaDir.empty() && opt(mPara, "conductivity")) mkdir_fs(sigmaDir);
    if (!AkwDir.empty() && opt(mPara, "Akw")) mkdir_fs(AkwDir);
    if (!pumpDir.empty() && opt(mPara, "pump")) mkdir_fs(pumpDir);
    // heisenberg dirs
    if (!SiSjDir.empty() && opt(mPara, "SiSj")) mkdir_fs(SiSjDir);
    if (!SkwDir.empty() && opt(mPara, "Skw")) mkdir_fs(SkwDir);
    if (!RamanDir.empty() && opt(mPara, "Raman")) mkdir_fs(RamanDir);
}

void Path::print(std::ostream& os, bool brief) const {
    printLine(50, '-');
    os << "root dir: " << rootDir << '\n';
    os << "project: " << project << '\n';
    os << "basis dir: " << basisDir << '\n';
    os << "job data dir :" << parameterDir << '\n';
    if (!brief) {
        printLine(50, '-');
        if (!wavefuncDir.empty()) {
            os << "wavefunc dir: " << wavefuncDir << '\n'; 
        }
        // hubbard files
        if (!particleDistrDir.empty()) {
            os << "particle distribution dir: " << particleDistrDir << '\n';
        }
        if (!sigmaDir.empty()) {
            os << "conductivity dir: " << sigmaDir << '\n'; 
        }
        if (!AkwDir.empty()) {
            os << "akw dir: " << AkwDir << '\n';
        }
        if (!pumpDir.empty()) {
            os << "pump dir: " << pumpDir << '\n'; 
        }
        // heisenberg files
        if (!SiSjDir.empty()) {
            os << "sisj dir: " << SiSjDir << '\n'; 
        }
        if (!SkwDir.empty()) {
            os << "skw dir: " << SkwDir << '\n';
        }
        if (!RamanDir.empty()) {
            os << "raman dir: " << RamanDir << '\n';
        }
    }
}

#endif // __PATH_H__