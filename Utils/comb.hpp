#pragma once

#include <iostream>
#include <cmath>

/*
    ***************
    * Combinatory *
    ***************
*/
unsigned long long pow(unsigned long long base, unsigned long long power);

long long pow(long long base, long long power);

template <class Type>
Type factorial(Type n) {
    if (n < 0) {
        std::cerr<<"n must be a non-negative integer!"; 
        exit(EXIT_FAILURE);
    }
    Type f = 1;
    while (n > 0) {
        f *= n;
        --n;
    }
    return f;
}

template <class Type>
Type combination(Type n, Type k) {
    if (n < k or n < 1) {
        std::cerr<<"n can't be smaller than k and both n,k should be positive integer!\n"; 
        exit(EXIT_FAILURE);
    }
    if (n == k or k == 0) return 1;
    Type n1 = k < (n-k) ? k : n-k;
    Type* v1 = new Type[n1];
    Type* v2 = new Type[n1];
    for (Type i = 0; i < n1; ++i) {
        v1[i] = i + 1;
        v2[i] = n - n1 + i + 1;
    }
    for (Type i = n1; i > 0; --i) {
        for (Type j = 0; j < n1; ++j) {
            if (v2[j] % v1[i-1] == 0) {
                v2[j] /= v1[i-1];
                v1[i-1] = 1;
                break;
            }
        }
    }
    Type result1 = 1;
    Type result2 = 1;
    for (Type i = 0; i < n1; ++i) {
        result1 *= v1[i];
        result2 *= v2[i];
    }
    delete[] v1;
    delete[] v2;
    return result2/result1;
}