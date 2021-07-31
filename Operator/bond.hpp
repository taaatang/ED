#ifndef __BOND_H__
#define __BOND_H__
/**
 * @file bond.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief describe lattice site ids involved in an interaction term.
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include "Global/globalType.hpp"
#include "Geometry/Transform.hpp"

/**
 * @brief sited ids form a bond
 * @example J*Si*Sj => bond.val = J, bond.sites = {i, j}
 * @tparam T bond val type
 * @tparam N bond size
 */
template<typename T, size_t N>
struct Bond {
	T val;
	std::array<int, N> sites;
public:
	Bond() = default;
	Bond(T val_, std::array<int, N> sites_) : val(val_), sites(sites_){ }
	inline int& operator[](int idx) { return sites.at(idx); }
	inline const int& operator[](int idx) const { return sites.at(idx); }
};

template<typename T, size_t N>
Bond<T,N> operator*(const Transform<T>& g, const Bond<T,N>& b) {
	Bond<T,N> gb;
	gb.val = b.val;
	for (size_t i = 0; i < N; ++i) {
		gb[i] = g[b[i]];
	}
	return gb;
}

template<typename T, size_t N>
bool operator==(const Bond<T, N>& lhs, const Bond<T, N>& rhs) {
	for (size_t i = 0; i < N; ++i) {
		if (lhs[i] != rhs[i]) {
			return false;
		}
	}
	return true;
}

template<typename T, size_t N>
bool operator<(const Bond<T, N>& lhs, const Bond<T, N>& rhs) {
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
	size_t i = 0;
	for (size_t j = 1; j < N; ++j) {
		if (b[j] < b[i]) {
			i = j;
		}
	}
	std::array<int, N> sorted;
	for (size_t j = 0; j < N; ++j) {
		sorted[j] = b[(i + j) % N];
	}
	b.sites = sorted;
}

#endif // __BOND_H__