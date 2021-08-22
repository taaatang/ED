#include <iomanip>
#include "utils/io.hpp"

std::string tostr(double val, int digit){
    std::ostringstream strTmp;
    strTmp<<std::fixed<<std::setprecision(digit)<<val;
    return strTmp.str();
}

std::string tostr(int val){
    return std::to_string(val);
}

void tolower(std::string &str) {
    for (auto &c : str) {
        c = tolower(c);
    }
}

void toupper(std::string &str) {
    for (auto &c : str) {
        c = toupper(c);
    }
}

void printLine(int n, char c) {
    std::cout << std::string(n, c) << '\n';
}