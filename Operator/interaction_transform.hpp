#ifndef __INTERACTION_TRANSFORM_H__
#define __INTERACTION_TRANSFORM_H__
/**
 * @file interaction_transform.hpp
 * @author Ta Tang (tatang.physics@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "Geometry/Generator.hpp"
#include "Operator/interaction.hpp"
#include "Operator/links.hpp"

/**
 * @brief g.factor * Op * g, where g is a symmetry operation.
 * 
 * @tparam T 
 * @tparam N 
 * @details When acting on a basis state |r>, first transform it to |gr>,  then apply Op: Op|gr>, 
 * where Op are interaction terms of size N
 */
template <typename T, size_t N>
struct TrInteractions {
	Transform<T> g;
	Interactions<T, N> Op;
	TrInteractions(const Transform<T>& g_, LINK_TYPE type):g(g_), Op(type) {}
};
/**
 * @brief TrInteractions are considered "equal" if they are associated with the same symmetry transformation
 * 
 * @tparam T 
 * @tparam N 
 * @param lhs 
 * @param rhs 
 * @return true 
 * @return false 
 */
template <typename T, size_t N>
inline bool operator==(const TrInteractions<T, N>& lhs, const TrInteractions<T, N>& rhs) {
	return lhs.g == rhs.g;
}

template <typename T, size_t N>
inline bool operator<(const TrInteractions<T, N>& lhs, const TrInteractions<T, N>& rhs) {
	return lhs.g < rhs.g;
}

template <typename T, size_t N>
Interactions<T, N> linkToTrInteractions(const Transform<T>& g, const Link<T>& link) {
	return g * linkToInteraction<T, N>(link);
} 

#endif // __INTERACTION_TRANSFORM_H__