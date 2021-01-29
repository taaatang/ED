#include "../Global/global.hpp"
#include "../Utils/paras.hpp"
#include "../Geometry/Geometry.hpp"
#include <memory>

std::string configFile = "../Input/config.txt";
Parameters para(configFile);
std::unique_ptr<Geometry> latt; // to be deleted ...
std::unique_ptr<Basis> Bi;
std::unique_ptr<Basis> Bf;
int nu;
int nd;

setlatt(para, latt);