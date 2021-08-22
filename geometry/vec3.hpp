#pragma once

#include <array>
#include <iostream>
#include <cmath>

template <typename T>
using Vec3 = std::array<T, 3>;

using Vec3d = Vec3<double>;
using Vec3i = Vec3<int>;

template <typename T>
inline Vec3<T> operator+(const Vec3<T>& v1, const Vec3<T>& v2) {
    Vec3<T> res;
    for (int i = 0; i < 3; ++i) {
        res[i] = v1[i] + v2[i];
    }
    return res;
}

template <typename T>
inline Vec3<T> operator-(const Vec3<T>& v1, const Vec3<T>& v2) {
    Vec3<T> res;
    for (int i = 0; i < 3; ++i) {
        res[i] = v1[i] - v2[i];
    }
    return res;
}

template <typename T>
inline Vec3<T> operator*(const Vec3<T>& v1, const Vec3<T>& v2) {
    Vec3<T> res;
    for (int i = 0; i < 3; ++i) {
        res[i] = v1[i] * v2[i];
    }
    return res;
}

template <typename T>
inline Vec3<T> operator/(const Vec3<T>& v1, const Vec3<T>& v2) {
    Vec3<T> res;
    for (int i = 0; i < 3; ++i) {
        res[i] = v1[i] / v2[i];
    }
    return res;
}

template <typename T>
inline bool operator==(const Vec3<T>& v1, const Vec3<T>& v2) {
    return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
}

template <typename T>
inline T dot(const Vec3<T>& v1, const Vec3<T>& v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template <typename T>
inline void normalize(Vec3<T>& v) {
    auto norm = std::sqrt(dot(v, v));
    v = v / norm;
}

template <typename T>
inline Vec3<T> operator*(T a, const Vec3<T>& v) {
    Vec3<T> res;
    for (int i = 0; i < 3; ++i) {
        res[i] = a * v[i];
    }
    return res;
}

template <typename T>
inline Vec3<T> operator*(const Vec3<T>& v, T a) {
    return a * v;
}

template <typename T>
inline Vec3<T> operator/(const Vec3<T>& v, T a) {
    return 1.0/a * v;
}

template <typename T>
inline std::ostream& operator<<(std::ostream &os, Vec3<T> v) {
    os<<"["<<v[0]<<" "<<v[1]<<" "<<v[2]<<"]";
    return os;
}