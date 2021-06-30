#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include <iostream>
#include <vector>
#include "Utils/io.hpp"
#include "Utils/bitop.hpp"

template <typename T>
struct Transform {
	T factor;
	size_t size;
	std::vector<int> transformList;
	Transform(T fac, size_t N):factor(fac), size(N), transformList(N, 0) { }
	//? why this caused trouble in copying Transform
	// Transform(const Transform<T>& t) { factor = t.factor; transformList = t.transformList;}
	int& operator[](int idx) { return transformList.at(idx); }
	const int& operator[](int idx) const {return transformList.at(idx); }
	idx_t tr(idx_t repI) const;
	idx_t tr(idx_t repI, double& sgn) const;
};

template<typename T>
idx_t Transform<T>::tr(idx_t repI) const {
	idx_t repF{0};
	for (int i = 0; i < (int)size; ++i) {
		if (bitTest(repI, i)) {
			bitSet(repF, transformList[i]);
		}
	}
	return repF;
}

template<typename T>
idx_t Transform<T>::tr(idx_t repI, double& sgn) const {
	idx_t repF{0};
	std::vector<int> seq;
	for (int i = 0; i < size; ++i) {
		if (bitTest(repI, i)) {
			bitSet(repF, transformList[i]);
			seq.push_back(transformList[i]);
		}
	}
	sgn = seqSign(seq);
	return repF;
}


template <typename T>
inline void check(const Transform<T>& lhs, const Transform<T>& rhs) {
	assert_msg(lhs.size == rhs.size, "mismatch transformation size:" + tostr(int(lhs.size)) + " vs " + tostr(int(rhs.size)));
}

template <typename T>
inline Transform<T> operator*(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	Transform<T> res(lhs.factor * rhs.factor, lhs.size);
	for (size_t i = 0; i < res.size; ++i) {
		res[i] = lhs[rhs[i]];
	}
	return res;
}

template <typename T>
bool operator==(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	for (size_t i = 0; i < lhs.size; ++i) {
		if (lhs[i] != rhs[i]) return false;
	}
	return true;
}

template <typename T>
bool operator<(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	for (int i = 0; i < (int)lhs.size; ++i) {
		if (lhs[i] < rhs[i]) {
			return true;
		} else if (lhs[i] > rhs[i]) {
			return false;
		}
	}
	return false;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Transform<T>& t) {
	os << "Transformation Character:" << t.factor << std::endl;
	for (size_t i = 0; i < t.size; ++i) {
		os << i << "->" << t[i] << '\t';
	}
	os << std::endl;
	return os;
}

#endif // __TRANSFORM_H__