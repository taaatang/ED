#ifndef __INTERACTION_TRANSFORM_H__
#define __INTERACTION_TRANSFORM_H__

#include "Geometry/Generator.hpp"
#include "Operator/interaction.hpp"
#include "Operator/links.hpp"

template <typename T, size_t N>
struct TrInteractions {
	Transform<T> g;
	Interactions<T, N> Op;
	TrInteractions(const Transform<T>& g_, LINK_TYPE type):g(g_), Op(type) {}
};

template <typename T, size_t N>
inline bool operator==(const TrInteractions<T, N>& lhs, const TrInteractions<T, N>& rhs) {
	return lhs.g == rhs.g;
}

template <typename T, size_t N>
inline bool operator<(const TrInteractions<T, N>& lhs, const TrInteractions<T, N>& rhs) {
	return lhs.g < rhs.g;
}

template <typename T, size_t N>
Interactions<T, N> trLink(const Transform<T>& g, const Link<T>& link) {
	assert_msg(link.getLinkSize() + 1 == N, "link size mismatch interaction bond size!");
	Interactions<T, N> Op(link.getLinkType());
	T val = link.getVal();
	for (const auto& sites : link.bond()) {
		Bond<T, N> b;
		b.val = val;
		for (size_t i = 0; i < N; ++i) {
			b[i] = g[sites.at(i)];
		}
		Op.add(b);
	}
	return Op;
} 


#endif // __INTERACTION_TRANSFORM_H__