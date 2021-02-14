#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cmath>
#include <memory>
#include "Global/globalPara.hpp"
#include "Utils/paras.hpp"
#include "Utils/path.hpp"
#include "Geometry/Geometry.hpp"
#include "Operator/OperatorBasics.hpp"
#include "Pulse/pulse.hpp"

// std::string configFile = "../Input/config.txt";
// Parameters para(configFile);
std::string lattstr, occstr, tpstr;  
std::string inputDir = "Input";
std::string rootDir, projectDir;
std::string basisDir, wavefuncDir, corrDir, SkwDir, AkwDir, sigmaDir, pumpDir;

Parameters pathPara(inputDir, {"path.txt"});
Parameters modelPara(inputDir, {"lattice.txt", "Hubbard.txt"});
Parameters pulsePara(inputDir, {"pulse.txt"});

Path path(&pathPara, &modelPara, &pulsePara);

std::unique_ptr<Geometry> latt;

std::unique_ptr<Basis> Bi, Bf;

std::unique_ptr<OperatorBase<dataType>> H, Hf;

Pulse pulse;


#endif // __CONFIG_H__