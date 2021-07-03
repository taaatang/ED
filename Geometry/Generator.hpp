#ifndef __GENERATOR_H__
#define __GENERATOR_H__
/**
 * @file Generator.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <algorithm>
#include <cmath>
#include "Geometry/Transform.hpp"

/**
 * 
 * @brief Generator for creating eigen states of a symmetry group G based on its 1-dimension irreducible representations.
 * 
 * @tparam T type of the character chi
 * 
 * @details  A particular representation repi = {(chi_g, g), g \in G},
 * Generator P_repi = 1/|G| * \sum_g chi_g * Ug, Ug is the corresponding symmetry transformation operator
 * Then acting on a basis P_repi|b> create eigen states for all Ug with eigen values conj(chi_g)
 */
template <typename T>
struct Generator {
	// symmetry group transformations and corresponding characters
	std::vector<Transform<T>> U{};
	// symmetry group size
	size_t G{0};
	bool condensed{false};
public:
	void setIdentity(int N);
	Transform<T>& operator[](int idx) { return U.at(idx); }
	const Transform<T>& operator[](int idx) const { return U.at(idx); }
	void add(const Transform<T>& u) { G++; U.push_back(u); }
	void condense();
};

// set the generator to identity 
template<typename T>
void Generator<T>::setIdentity(int N) {
	U.clear();
	G = 0;
	Transform<T> t(1.0, N);
	for (int i = 0; i < N; ++i) {
		t[i] = i;
	}
	this->add(t);
}

/**
 * @brief Combine same transformations: if Ug' = Ug, then set chi_g += chi_g' and remove Ug'
 * 
 * @tparam T 
 */
template <typename T>
void Generator<T>::condense() {
	std::sort(U.begin(), U.end());
	auto cur = U.begin();
	auto next = cur + 1;
	while (next != U.end()) {
		if (*next == *cur) {
			cur->factor += next->factor;
		} else {
			cur = next;
		}
		next++;
	}
	auto end = std::unique(U.begin(), U.end());
	if (end != U.end()) {
		// std::cout << "original size:" << G << std::endl;
		U.erase(end, U.end());
		G = U.size();
		// std::cout << "condensed size:" << G << std::endl;
	}
	condensed = true;
}

template <typename T>
Generator<T> operator*(const Generator<T>& lhs, const Generator<T>& rhs) {
	Generator<T> res;
	for (size_t i = 0; i < lhs.G; ++i) {
		for (size_t j = 0; j < rhs.G; ++j) {
			res.add(lhs[i] * rhs[j]);
		}
	}
	res.condense();
	return res;
}

/**
 * @brief equal only if for each symmetry element g: chi_g and Ug are equal
 * 
 * @tparam T 
 * @param lhs 
 * @param rhs 
 * @return true 
 * @return false 
 */
template <typename T>
bool operator==(const Generator<T>& lhs, const Generator<T>& rhs) {
	assert_msg(lhs.condensed && rhs.condensed, "Generators should be condensed first before comparison!");
	if (lhs.G != rhs.G) return false;
	for (int i = 0; i < lhs.G; ++i) {
		if (lhs[i] == rhs[i] and std::abs(lhs[i].factor - rhs[i].factor) < 1e-8) {
			continue;
		} else {
			return false;
		}
	}
	return true;
}

template <typename T>
bool commute(const Generator<T>& lhs, const Generator<T>& rhs) {
	auto lr = lhs * rhs;
	auto rl = rhs * lhs;
	return lr == rl;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Generator<T>& g) {
	os << "Generator Size:" << g.G << std::endl;
	for (size_t i = 0; i < g.G; ++i) {
		os << g[i];
	}
	return os;
}

#endif // __GENERATOR_H__