#include <iostream>
#include "Solver/BiCGSTAB.hpp"
#include "Measure/measure.hpp"

int main( ) {
	int workerID, workerNum;
	bool isMaster;
	init(workerID, workerNum, isMaster);
	int dim = 100;
	Identity<cdouble> eye(dim);
	cdouble z{2.0, 0.5};
	std::vector<cdouble> x(dim, 0.0);
	std::vector<cdouble> b(dim, 3.0);
	int iterCount{0};
	double res{0.0};
	// BiCGSTAB(&eye, z, b.data(), x.data(), eye.nloc, &eye, iterCount, res);
	// std::cout << "iter: " << iterCount << ", err: " << res << '\n' << "solution:\n" << x << '\n'; 

	Basis B(dim);
	SparseMatrix<cdouble> A(&B, &B, 1, 0);
	cdouble a = 1.0;
	cdouble c = 2.0;
	MAP<cdouble> frow({{idx_t(0), a}, {idx_t(1), c}});
	A.pushRow(&frow);
	for (idx_t i = 1; i < dim -1; ++i) {
		MAP<cdouble> row({{i - 1, c}, {i, a}, {i + 1, c}});
		A.pushRow(&row);
	}
	MAP<cdouble> lrow({{idx_t(dim - 2), c}, {idx_t(dim - 1), a}});
	A.pushRow(&lrow);
	A.build();
	scale(x.data(), cdouble(0.0), dim);
	BiCGSTAB(&A, z, b.data(), x.data(), A.nloc, &eye, iterCount, res);
	printLine();
	std::cout << "iter: " << iterCount << ", err: " << res << '\n' << "solution:\n" << x << '\n';
	std::vector<cdouble> bp(dim, 0.0);
	A.MxV(x.data(), bp.data());
	axpy(bp.data(), z, x.data(), dim);
	printLine();
	std::cout << "b:\n" << b << '\n';
	printLine();
	std::cout<< "Ax:\n" << bp << '\n';	
	MPI_Finalize();
    return 0;
}