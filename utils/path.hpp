#pragma once

#include "utils/paras.hpp"
// #include "utils/io.hpp"

class Path {
public:
    Path() = default;
    Path(const Parameters* path, const Parameters* model, const Parameters* pulse = nullptr);
    
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

    std::string NeqNelDir;
    std::string NeqNphDir;
    std::string NeqSkwDir;
};