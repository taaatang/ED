#ifndef __LINKS_H__
#define __LINKS_H__

#include <vector>
#include "Operator/OperatorBasics.hpp"
#include "Global/globalPara.hpp"

std::vector<Link<dataType>> HeisenbergLink( );

std::vector<Link<dataType>> HubbardSingleBandLink( );

std::vector<Link<dataType>> HubbardMultiBandLink( );

std::vector<Link<dataType>> HubbardLink( );

#endif // __LINKS_H__