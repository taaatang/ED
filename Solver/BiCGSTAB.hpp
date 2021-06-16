#ifndef __BICGSTAB_H__
#define __BICGSTAB_H__

#include <vector>
#include <cmath>

#include "Global/globalType.hpp"
#include "Operator/SparseMatrix.hpp"
#include "Algebra/algebra.hpp"

// solve (A + z) x = b left precomditioned by  Minv * (A + z) ~ I. Minv is the simple Jacobi preconditioner.
// based on Numerical Linear Algebra with Julia, Eric Darve, chp9.4; and wikipedia BiCGSTAB page.
template <class T>
void BiCGSTAB(BaseMatrix<T> *A, T z, const T *b, T *x, idx_t size, BaseMatrix<T> *Minv, int &iterCount, double &res, int iterMax = 500, double tol = 1e-8, double zero = 1e-15) {
	res = 0.0;
	auto nb = mpiNorm(b, size);
	// r0 = b - Ax
	auto r0 = std::vector<T>(size, 0.0);
	auto s = std::vector<T>(size, 0.0);
	A->MxV(x, r0.data());
	axpy(r0.data(), z, x, size);
	scale(r0.data(), -T(1.0), size);
	axpy(r0.data(), T(1.0), b, size);
	// rh0 = r0
	auto rh0 = std::vector<T>(size, 0.0);
	copy(size, r0.data(), rh0.data());

	T rho_0{1.0}, alpha{1.0}, w0{1.0};

	auto v0 = std::vector<T>(size, 0.0);
	auto p0 = std::vector<T>(size, 0.0);

	int restart = 0;
	for (iterCount = 0; iterCount < iterMax; ++iterCount) {
		// BiCG
		auto rho_1 = mpiDot(rh0.data(), r0.data(), size);
		if (std::abs(rho_0) < zero || std::abs(w0) < zero) {
			// std::cout << "BiCG stoped at step " << iterCount << '\n' << "rho: " << rho_0 << ", w: " << w0 << ", res: " << res << '\n';
			// exit(1);
			copy(size, r0.data(), rh0.data(), size);
			rho_1 = mpiDot(rh0.data(), r0.data(), size);
			restart++;
			if (restart == 1) iterCount = 0;	
			return;
		}
		auto beta = (rho_1 / rho_0) * (alpha / w0);
		// pi = ri-1 + beta * (pi-1 - wi-1 * vi-1)
		axpy(p0.data(), -w0, v0.data(), size);
		combine(p0.data(), T(1.0), r0.data(), beta, p0.data(), size);
		// v0 = A * K2 * K1 * p0
		auto y = std::vector<T>(size, 0.0);
		Minv->MxV(p0.data(), y.data());
		A->MxV(y.data(), v0.data());
		axpy(v0.data(), z, y.data(), size);
		alpha = rho_1 / mpiDot(rh0.data(), v0.data(), size);
		axpy(x, alpha, y.data(), size);
		combine(s.data(), T(1.0), r0.data(), -alpha, v0.data(), size);
		res = mpiNorm(s.data(), size) / nb;
		// resVec.push_back(res);
		if (res < tol) {
			return;
		}

		// GMRES
		auto sh = std::vector<T>(size, 0.0);
		Minv->MxV(s.data(), sh.data());
		auto t = std::vector<T>(size, 0.0);
		A->MxV(sh.data(), t.data());
		axpy(t.data(), z, sh.data(), size);
		w0 = mpiDot(t.data(), s.data(), size) / mpiDot(t.data(), t.data(), size);
		axpy(x, w0, sh.data(), size);
		combine(r0.data(), T(1.0), s.data(), -w0, t.data(), size);
		res = mpiNorm(r0.data(), size) / nb;
		// resVec.push_back(res);
		if (res < tol) {
			return;
		}
		rho_0 = rho_1;
	}	
} 
#endif // __BICGSTAB_H__