#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <algorithm>
#include <cmath>
#include "Geometry/Transform.hpp"

template <typename T>
struct Generator {
	std::vector<Transform<T>> U{};
	size_t G{0};
	bool condensed{false};
	Generator() { }
	void identity(int N);
	Transform<T>& operator[](int idx) { return U.at(idx); }
	const Transform<T>& operator[](int idx) const { return U.at(idx); }
	void add(const Transform<T>& u) { G++; U.push_back(u); }
	void condense();
};

template<typename T>
void Generator<T>::identity(int N) {
	U.clear();
	G = 0;
	Transform<T> t(1.0, N);
	for (int i = 0; i < N; ++i) {
		t[i] = i;
	}
	this->add(t);
}


template <typename T>
void Generator<T>::condense() {
	if (!condensed) {
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

template <typename T>
bool operator==(Generator<T>& lhs, Generator<T>& rhs) {
	lhs.condense();
	rhs.condense();
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
bool commute(Generator<T>& lhs, Generator<T>& rhs) {
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