#ifndef __INTERACTION_H__
#define __INTERACTION_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include "Operator/bond.hpp"
#include "Geometry/Generator.hpp"

template<typename T, size_t N>
struct Interactions {
	LINK_TYPE type;
	std::vector<Bond<T, N>> bonds{};
	Interactions(LINK_TYPE t) : type(t) { }
	void add(const Bond<T, N>& b) { bonds.push_back(b); }
	void condense(char order = 'n');
	Interactions<T, N>& operator+=(const Interactions<T, N>& rhs) {
		assert_msg(type == rhs.type, "Interactions type mismatch in operator+=!");
		bonds.insert(bonds.end(), rhs.bonds.cbegin(), rhs.bonds.cend());
		return *this;
	}
};

template<typename T, size_t N>
Interactions<T, N> operator*(const Interactions<T, N>& op, T fac) {
	Interactions<T, N> res(op);
	for (auto& bond : res.bonds) {
		bond.val *= fac;
	}
	return res;
}

template<typename T, size_t N>
Interactions<T, N> operator*(T fac, const Interactions<T, N>& op) {
	return op * fac;
}

template<typename T, size_t N>
void Interactions<T, N>::condense(char order) {
	if (order == 'n') {
		for (auto& b : bonds) {
			normalOrder(b);
		}
	} else if (order == 'c') {
		for (auto& b : bonds) {
			cyclicOrder(b);
		}
	}
	std::sort(bonds.begin(), bonds.end());
	auto cur = bonds.begin();
	auto next = cur + 1;
	while (next != bonds.end()) {
		if (*next == *cur) {
			cur->val += next->val;
		} else {
			cur = next;
		}
		next++;
	}
	auto end = std::unique(bonds.begin(), bonds.end());
	if (end != bonds.end()) {
		bonds.erase(end, bonds.end());
	}
	for (auto last = bonds.end() - 1; last >= bonds.begin(); last--) {
		if (std::abs(last->val) < 1e-8) {
			last = bonds.erase(last);
		}
	}
}

// assumming lhs and rhs are condensed
template<typename T, size_t N>
bool operator==(const Interactions<T, N>& lhs, const Interactions<T, N>& rhs) {
	if (lhs.type != rhs.type) {
		return false;
	}
	if (lhs.bonds.size() != rhs.bonds.size()) {
		return false;
	}
	for (int i = 0; i < lhs.bonds.size(); ++i) {
		if (!(lhs.bonds[i] == rhs.bonds[i]) or (std::abs(lhs.bonds[i].val - rhs.bonds[i].val) > 1e-8)) {
			return false;
		}
	}
	return true;
}


#endif // __INTERACTION_H__