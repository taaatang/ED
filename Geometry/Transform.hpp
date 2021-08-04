#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include <iostream>
#include <vector>
#include "Utils/io.hpp"
#include "Utils/bitop.hpp"

/**
 * @brief lattice symmetry transformation factor * g, where g is a symmetry operation (translation or point group symm). 
 * g takes orbital i to transformList[i] 
 * 
 * @tparam T type of factor; should generally be complex
 */
template <typename T>
struct Transform {
	T factor;
	size_t size;
	std::vector<int> transformList;
public:
	Transform(T fac, size_t N):factor(fac), size(N), transformList(N, 0) { }
	int& operator[](int idx) { return transformList.at(idx); }
	const int& operator[](int idx) const {return transformList.at(idx); }
	idx_t tr(idx_t repI) const;
	idx_t tr(idx_t repI, double& sgn) const;
	std::vector<uint8_t> tr(const std::vector<uint8_t>& rep) const;
};
/**
 * @brief transform basis sate repI to g(repI)
 * 
 */
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
/**
 * @brief transform basis repI to g(repI) and returns the fermion permutation sign
 * 
 * @tparam T 
 * @param repI input basis state
 * @param sgn fermion sign due to the permutation of the transformation
 * @return idx_t transformed basis state 
 */
template<typename T>
idx_t Transform<T>::tr(idx_t repI, double& sgn) const {
	idx_t repF{0};
	std::vector<int> seq;
	for (size_t i = 0; i < size; ++i) {
		if (bitTest(repI, i)) {
			bitSet(repF, transformList[i]);
			seq.push_back(transformList[i]);
		}
	}
	sgn = seqSign(seq);
	return repF;
}

template<typename T>
std::vector<uint8_t> Transform<T>::tr(const std::vector<uint8_t>& rep) const {
	std::vector<uint8_t> res(rep.size(), 0);
	for (size_t i = 0; i < size; ++i) {
		if (rep.at(i)) {
			res.at(transformList[i]) = rep[i];
		}
	}
	return res;
}


template <typename T>
inline void check(const Transform<T>& lhs, const Transform<T>& rhs) {
	assert_msg(lhs.size == rhs.size, "mismatch transformation size:" + tostr(int(lhs.size)) + " vs " + tostr(int(rhs.size)));
}
/**
 * @brief combines two transformations g_left*g_right
 * 
 * @tparam T 
 * @param lhs 
 * @param rhs 
 * @return Transform<T> 
 */
template <typename T>
inline Transform<T> operator*(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	Transform<T> res(lhs.factor * rhs.factor, lhs.size);
	for (size_t i = 0; i < res.size; ++i) {
		res[i] = lhs[rhs[i]];
	}
	return res;
}
/**
 * @brief two transformations are equal if their transfomationLists are the same
 * 
 */
template <typename T>
bool operator==(const Transform<T>& lhs, const Transform<T>& rhs) {
	check(lhs, rhs);
	for (size_t i = 0; i < lhs.size; ++i) {
		if (lhs[i] != rhs[i]) return false;
	}
	return true;
}
/**
 * @brief compare the transformationList from the first entry to the last entry
 * 
 */
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