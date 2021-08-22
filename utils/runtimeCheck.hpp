#pragma once

#include <cassert>
#include <string>
#include <iostream>

#define LOCATION(cond) if (cond) {std::cout << "line number: "<< __LINE__ << "\nfrom func: " << __func__ << "\nfrom file: "<< __FILE__ << '\n';}

inline void exit_msg(const std::string& msg){
    std::cerr << msg <<std::endl;
    exit(EXIT_FAILURE);
}

inline void assert_msg(bool condition, const std::string& msg){
    if (!condition) {
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }
}