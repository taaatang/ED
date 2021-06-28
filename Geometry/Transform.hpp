#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include <iostream>
#include <vector>
#include "Utils/io.hpp"

template <typename T>
struct Transform {
	T factor;
	size_t size;
	std::vector<int> transformList;
	Transform(T fac, size_t N):factor(fac), size(N), transformList(N, 0) { }
	int& operator[](int idx) { return transformList.at(idx); }
};

template <typename T>
inline void check(const Transform<T>& lhs, const Transform<T>& rhs) {
	assert_msg(lhs.size == rhs.size, "mismatch transformation size!");
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
	for (int i = 0; i < lhs,size; ++i) {
		if (lhs[i] != rhs[i]) return false;
	}
	return true;
}

template <typename T>
bool operator<(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	for (int i = 0; i < lhs.size; ++i) {
		if (lhs[i] < rhs[i]) return true;
		if (lhs[i] > rhs[i]) return false;
	}
	return false;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Transform<T>& t) {
	os << "Transformation Character:" << t.factor << std::endl;
	for (int i = 0; i < N; ++i) {
		os << i << "->" << t[i] << '\t';
	}
	os << std::endl;
	return os;
}

#endif // __TRANSFORM_H__