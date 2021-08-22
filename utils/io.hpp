/**
 * @file io.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-08-22
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>

inline std::string tostr(double val, int digit = 2){
    std::ostringstream strTmp;
    strTmp << std::fixed << std::setprecision(digit) << val;
    return strTmp.str();
}

inline std::string tostr(int val){
    return std::to_string(val);
}

inline void tolower(std::string &str) {
    for (auto &c : str) {
        c = tolower(c);
    }
}

inline void toupper(std::string &str) {
    for (auto &c : str) {
        c = toupper(c);
    }
}

inline void printLine(int n = 50, char c = '*') {
    std::cout << std::string(n, c) << '\n';
}

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