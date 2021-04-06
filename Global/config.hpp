#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cmath>
#include <memory>
#include "Global/globalPara.hpp"
#include "Utils/paras.hpp"
#include "Utils/path.hpp"
#include "Geometry/Geometry.hpp"

// std::string configFile = "../Input/config.txt";
// Parameters para(configFile);
std::string inputDir = "heisInput";

Parameters pathPara(inputDir, {"path.txt"});
Parameters modelPara(inputDir, {"lattice.txt", "hamiltonian.txt"});
Parameters pulsePara(inputDir, {"pulse.txt"});
Parameters measurePara(inputDir, {"measure.txt"});

Path path(&pathPara, &modelPara, &pulsePara);

std::unique_ptr<Geometry> latt;

std::unique_ptr<Basis> Bi, Bf;

std::unique_ptr<HamiltonianBase<dataType>> H, Hf;


#endif // __CONFIG_H__