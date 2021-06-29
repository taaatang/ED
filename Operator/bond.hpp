#ifndef __BOND_H__
#define __BOND_H__

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "Global/globalType.hpp"

template<typename T, size_t N>
struct Bond {
	T val;
	std::array<int, N> sites;
	inline int& operator[](int idx) { reurn sites.at(idx); }
	inline const int& operator[](int idx) const { return sites.at(idx); }
};

template<typename T, size_t N>
bool operator==(const Bond<T, N>& lhs, const Bond<T, N> rhs) {
	for (size_t i = 0; i < N; ++i) {
		if (lhs[i] != rhs[i]) {
			return false;
		}
	}
	return true;
}

template<typename T, size_t N>
bool operator<(const Bond<T, N>& lhs, const Bond<T, N> rhs) {
	for (size_t i = 0; i < N; ++i) {
		if (lhs[i] < rhs[i]) {
			return true;
		} else if (lhs[i] > rhs[i]) {
			return false;
		}
	}
	return false;
}

template<typename T, size_t N>
void normalOrder(Bond<T, N>& b) {
	std::sort(b.sites.begin(), b.sites.end());
}

template<typename T, size_t N>
void cyclicOrder(Bond<T, N>& b) {
	int i = 0;
	for (int j = 1; j < N; ++j) {
		if (b[j] < b[i]) {
			i = j;
		}
	}
	array<int, N> sorted;
	for (int j = 0; j < N; ++j) {
		sorted[j] = b[(i + j) % N];
	}
	b.sites = sorted;
}

template<typename T, size_t N>
struct Interactions {
	LINK_TYPE type;
	std::vector<Bond<T, N>> bonds;
}

#endif // __BOND_H__