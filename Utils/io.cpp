#include <iomanip>
#include "utils/io.hpp"
#include "global/constant.hpp"


/******
 * IO *
 ******/

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

std::ostream& operator<<(std::ostream& os, const VecD& vec) {
    os<<"[";
    for (auto val:vec) {
        os<<" "<<val;
    }
    os<<" ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<cdouble>& vec) {
    os<<"[";
    for (auto val:vec) {
        os << " " << (std::abs(std::real(val)) > 1e-12 ? std::real(val) : 0) << "+" << (std::abs(std::imag(val)) > 1e-12 ? std::imag(val) : 0) << "i";
    }
    os<<" ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const VecI& vec) {
    os<<"[";
    for (auto val:vec) {
        os<<" "<<val;
    }
    os<<" ]";
    return os;
}