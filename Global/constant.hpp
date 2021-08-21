#pragma once
//
// constant.hpp
// ED
//
// Created by tatang on 11/19/19
//
#include <complex>
#include <numbers>

inline constexpr int MPI_MASTER = 0;

inline constexpr double INFINITESIMAL = 1e-8;

inline constexpr double PI = std::numbers::pi_v<double>;

inline constexpr std::complex<double> CPLX_I(0.0, 1.0);