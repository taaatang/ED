#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cmath>
#include <memory>
#include "Global/globalPara.hpp"
#include "Utils/paras.hpp"
#include "Geometry/Geometry.hpp"
#include "Operator/OperatorBasics.hpp"

std::string configFile = "../Input/config.txt";
Parameters para(configFile);
std::unique_ptr<Geometry> latt;
std::unique_ptr<Basis> Bi;
std::unique_ptr<Basis> Bf;


#endif // __CONFIG_H__