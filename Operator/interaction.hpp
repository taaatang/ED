/**
 * @file interaction.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef __INTERACTION_H__
#define __INTERACTION_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include "Geometry/Generator.hpp"
#include "Operator/bond.hpp"
#include "Operator/links.hpp"

/**
 * @brief struct contains all bonds for an interaction
 * 
 * @tparam T 
 * @tparam N number of sites involved in the interaction
 */
template<typename T, size_t N>
struct Interactions {
	LINK_TYPE type;
	std::vector<Bond<T, N>> bonds{};
	// delete default constructor: must specify interaction type
	Interactions() = delete;
	Interactions(LINK_TYPE t) : type(t) { }
	void add(const Bond<T, N>& b) { bonds.push_back(b); }
	void clear() { bonds.clear(); }
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

template <typename T, size_t N>
Interactions<T, N> operator*(const Transform<T>& g, const Interactions<T, N>& Op) {
	Interactions<T, N> trOp(Op.type);
	for (const auto& b : Op.bonds) {
		trOp.add(g * b);
	}
	return trOp;
}

/**
 * @brief sort and combine duplicated interaction terms
 * 
 * @tparam T 
 * @tparam N 
 * @param order 'n':the order of sites in a bond is irrelavant, eg., 123 = 132 = 213 = 231 = 312 = 321, set all bonds to "normal" order (smallest)
 * 	123 in this example before sorting.
 * 	'c':only cyclic orders are considered equivalent, 123 = 231 = 312 != 132 ...
 * 	if not 'n' or 'c', different orderings are considered to be different interaction terms
 */
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

template<typename T, size_t N>
Interactions<T, N> linkToInteraction(const Link<T>& link) {
	assert_msg(link.getOrbNum() == N, "link size mismatch interaction bond size!");
	Interactions<T, N> Op(link.getLinkType());
	T val = link.getVal();
	for (const auto& sites : link.bond()) {
		Bond<T, N> b;
		b.val = val;
		for (size_t i = 0; i < N; ++i) {
			b[i] = sites.at(i);
		}
		Op.add(b);
	}
	return Op;
}

#endif // __INTERACTION_H__