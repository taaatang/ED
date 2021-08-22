#pragma once

#include <iostream>
#include <vector>
#include <complex>

std::string tostr(double val, int digit = 2);

std::string tostr(int val);

void tolower(std::string &str);

void toupper(std::string &str);

void printLine(int n = 50, char c = '*');

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os<<"[";
    for (const auto& val : vec) {
        os<<" "<<val;
    }
    os<<" ]";
    return os;
}

template<typename T>
inline std::ostream& operator<<(std::ostream& os, std::complex<T> val) {
	os << " " << (std::abs(std::real(val)) > 1e-12 ? std::real(val) : 0) << "+" << (std::abs(std::imag(val)) > 1e-12 ? std::imag(val) : 0) << "i";
	return os;
}