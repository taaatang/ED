#include "utils/path.hpp"
#include "utils/readWrite.hpp"
#include "utils/io.hpp"

Path::Path(const Parameters* path, const Parameters* model, const Parameters* pulse) {
    rootDir = *path->template get<Str>("rootDir");
    project = *path->template get<Str>("project");

    lattice = *model->template get<Str>("lattice type") + sep + *model->template get<Str>("boundary condition");
    std::string lattSize;
    int lx = *model->template get<int>("lx");
    int ly = *model->template get<int>("ly");
    lattSize = tostr(lx);
    if (ly > 0) lattSize += "x" + tostr(ly);
    lattice += sep + lattSize;

    particleNum = tostr(*model->template get<int>("nu")) + "u" + tostr(*model->template get<int>("nd")) + "d" + tostr(*model->template get<int>("maxPhPerSite")) + "ph";

    lattSymm = "k" + tostr(*model->template get<int>("kidx"));
    auto pidx = *model->template get<int>("pidx");
    if (pidx >= 0) {
        lattSymm += "p" + tostr(pidx);
    }

    symmLabel = particleNum + sep + lattSymm;

    parameterLabel = *model->template get<Str>("parameter label");
    parameterDir = rootDir + sep + project + sep + lattice + sep + symmLabel + sep + parameterLabel;
    parameterFile = parameterDir + sep + "parameter.txt";

    basisDir = rootDir + sep + project + sep + lattice + sep + "Basis" + sep + symmLabel;
    basisFile = basisDir + sep + "basis";
    basisNormFile = basisDir + sep + "norm";

    wavefuncDir = parameterDir + sep + "wavefunc";
    evalFile = wavefuncDir + sep + "eval";
 
    particleDistrDir = parameterDir + sep + "occpt";

    sigmaDir = parameterDir + sep + "sigma";

    AkwDir = parameterDir + sep + "Akw";

    BkwDir = parameterDir + sep + "Bkw";

    SiSjDir = parameterDir + sep + "sisj";
    SiSjFile = SiSjDir + sep + "sisj";
    SkwDir = parameterDir + sep + "skw";
    RamanDir = parameterDir + sep + "raman";

    if (pulse) {
        pumpDir = parameterDir + sep + "pump" + \
        sep + "pol_" + *pulse->template get<Str>("polarization label") + \
        sep + "i_" + tostr(*pulse->template get<double>("fluence")) + \
        sep + "w_" + tostr(*pulse->template get<double>("frequency")) + \
        sep + "phi_" + tostr(*pulse->template get<double>("phase"));

        pumpFile = pumpDir + sep + "pumpParameter.txt";
        NeqNelDir = pumpDir + sep + "nel";
        NeqNphDir = pumpDir + sep + "nph";
        NeqSkwDir = pumpDir + sep + "skw";
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
    if (!parameterDir.empty()) {
		mkdir_fs(parameterDir);
	}
    if (!basisDir.empty() && opt(mPara, "basis")) {
		mkdir_fs(basisDir);
	}
    if (!wavefuncDir.empty() && opt(mPara, "state")) {
		mkdir_fs(wavefuncDir);
	}
    if (!sigmaDir.empty() && opt(mPara, "occ")) {
		mkdir_fs(particleDistrDir);
	}
    if (!sigmaDir.empty() && opt(mPara, "conductivity")) {
		mkdir_fs(sigmaDir);
	}
    if (!AkwDir.empty() && opt(mPara, "Akw")) {
		mkdir_fs(AkwDir);
	}
    if (!BkwDir.empty() && opt(mPara, "Bkw")) {
		mkdir_fs(BkwDir);
	}
    // heisenberg dirs
    if (!SiSjDir.empty() && opt(mPara, "SiSj")) {
		mkdir_fs(SiSjDir);
	}
    if (!SkwDir.empty() && opt(mPara, "Skw")) {
		mkdir_fs(SkwDir);
	}
    if (!RamanDir.empty() && opt(mPara, "Raman")) {
		mkdir_fs(RamanDir);
	}
    if (!pumpDir.empty() && opt(mPara, "pump")) {
		mkdir_fs(pumpDir);
        if(!NeqNelDir.empty() && opt(mPara, "neqOcc")) {
            mkdir_fs(NeqNelDir);
        }
        if(!NeqNphDir.empty() && opt(mPara, "neqOcc")) {
            mkdir_fs(NeqNphDir);
        }
        if(!NeqSkwDir.empty() && opt(mPara, "neqSkw")) {
            mkdir_fs(NeqSkwDir);
        }
	}
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
        if (!BkwDir.empty()) {
            os << "bkw dir: " << BkwDir << '\n';
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