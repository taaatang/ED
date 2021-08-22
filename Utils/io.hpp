#pragma once

#include <iostream>

#include "global/typeAlias.hpp"

std::string tostr(double val, int digit = 2);

std::string tostr(int val);

void tolower(std::string &str);

void toupper(std::string &str);

void printLine(int n = 50, char c = '*');

std::ostream& operator<<(std::ostream& os, const VecD& vec);
std::ostream& operator<<(std::ostream& os, const VecI& vec);
std::ostream& operator<<(std::ostream& os, const std::vector<cdouble>& vec);