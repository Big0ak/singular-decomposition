#pragma once

#include <unordered_map>
#include "math/math_utils.h"
#include "math/solving_equation.h"
#include "math/matrix.h"

#include "eigenvalues/Jacobi.h"

enum Solving
{
	Exact_solution,
	Approximate_solution
};

class SVD
{
public:
	SVD();
	SVD(Solving solving);
	~SVD();

	bool solve(Matrix A, int n, int m, Matrix& U, Matrix& S, Matrix& V);

private:
	void SVD_2x2();
	void SVD_3x3();
	void SVD_4x4();
	void SVD_nxn();
	void SVD_matrix(Matrix A, int n, int m, std::vector<double> eigenvalues, Matrix eigenvectors, Matrix& U, Matrix& S, Matrix& V);
	bool Addition_base(Matrix& A, int n, int m, int r);
	void normalization(Matrix& eigenvectors, int n);
};