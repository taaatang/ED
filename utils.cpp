//
//  utils.cpp
//  ED
//
//  Created by tatang on 10/27/19.
//  Copyright © 2019 tatang. All rights reserved.
//

#include "utils.hpp"

void errExit(std::string msg){
    std::cerr<<msg<<std::endl;
    exit(EXIT_FAILURE);
}

unsigned long long pow(unsigned long long base, unsigned long long power){
    unsigned long long result = 1;
    for (unsigned long long i = 0; i < power; i++){
        result *= base;
    }
    return result;
}

long long pow(long long base, long long power){
    long long result = 1;
    for (long long i = 0; i < power; i++){
        result *= base;
    }
    return result;
}

