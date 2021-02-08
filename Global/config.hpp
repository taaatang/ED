#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cmath>
#include <memory>
#include "Global/globalPara.hpp"
#include "Utils/paras.hpp"
#include "Geometry/Geometry.hpp"
#include "Operator/OperatorBasics.hpp"
#include "Pulse/pulse.hpp"

// std::string configFile = "../Input/config.txt";
// Parameters para(configFile);
std::string inputDir = "Input";
Parameters para(inputDir, {"lattice.txt", "Hubbard.txt"});
Parameters pulsePara(inputDir, {"pulse.txt"});
std::unique_ptr<Geometry> latt;
std::unique_ptr<Basis> Bi, Bf;
std::unique_ptr<OperatorBase<dataType>> H;
Pulse pulse;


#endif // __CONFIG_H__