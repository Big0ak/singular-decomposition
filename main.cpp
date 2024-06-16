#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <complex>
#include <algorithm>
#include <cstring>

//#include "math/solving_equation.h"
#//include "tests/test.h"
//#include "math/matrix.h"
//#include "math/math_utils.h"
#include "svd.h"
#include "eigenvalues/Jacobi.h"

#define matrix_double std::vector<std::vector<double>>

// test svd
int main()
{
	int n = 0, m = 0;

	n = 3; m = 3;
	Matrix A({ {1,1,3}, {1,5,1}, {3,1,1} });

	Matrix U, S, V;
	SVD svd;
	svd.solve(A, n, m, U, S, V);

	std::cout << "input" << std::endl;
	std::cout << A << std::endl;

	std::cout << "SVD: " << std::endl;
	std::cout << "U: " << std::endl;
	std::cout << U << std::endl;
	std::cout << "S: " << std::endl;
	std::cout << S << std::endl;
	std::cout << "V: " << std::endl;
	std::cout << V << std::endl;

	Matrix CHECK;
	CHECK = V * S * U.transposition();

	std::cout << "CHECK" << std::endl;
	std::cout << CHECK << std::endl;

	return 0;
}

// test Jacobi
int main2()
{
	// -----------------------------------------------------------------------
	// Ferrari test
	//std::vector<std::complex<double>> answer;
	//std::vector<double> poly;
	//poly.resize(5);
	//poly[0] = -10;
	//poly[1] = 16;
	//poly[2] = -5;
	//poly[3] = 3;
	//poly[4] = 2;

	//Equation::Ferrari_real_number(poly, answer);
	// -----------------------------------------------------------------------

	Jacobi jacobi(Matrix(matrix_double{ { 1, 1, 3 }, { 1,5,1 }, { 3,1,1 } }), 3);
	double epsilon_elem = 0.0001;
	std::vector<double> eigenvalues_jacobi;
	Matrix T;

	jacobi.solve(epsilon_elem, 20);
	eigenvalues_jacobi = jacobi.get_eigenvalues();
	T = jacobi.get_eigenvectors();

	return 0;
}

// test solve equation
int main3()
{
	//int n = 4, m = 3;
	//std::vector<double> a {5, 0, 10, 6};
	//std::vector<double> b {1, 2, 4};
	//std::vector<double> r = multiply_polynomials(a, b, n, m);

	int cols = 3, rows = 3;
	std::vector<std::vector<double>> m;
	m.resize(cols);
	for (int i = 0; i < m.size(); i++)
		m[i].resize(rows);

	m = { {1,1,3}, {1,5,1}, {3,1,1} };
	//m = { { 4, 3, 2 }, { -2,5,0 }, { 2,0,6 } };

// -----------------------------------------------------------------------
	//std::vector<double> lambda_poly = multiply_polynomials(
	//	multiply_polynomials({ m[0][0], -1 }, 2, { m[1][1], -1 }, 2),
	//	3,
	//	{ m[2][2], -1 },
	//	2
	//);

	//lambda_poly[0] += m[1][0] * m[2][1] * m[0][2] + m[0][1] * m[1][2] * m[2][0];

	////lambda_poly[0] = 30;
	////lambda_poly[1] = -19;
	////lambda_poly[2] = 0;
	////lambda_poly[3] = 1;

	//std::vector<std::complex<double>> answer;
	//Cardano_real_number(lambda_poly, answer);
// -----------------------------------------------------------------------

	Matrix a_t(m);
	a_t.transposition();

	Matrix a(m);

	Matrix gram;
	gram = a_t * a;

	std::cout << gram;

	// Cardano
	/*
	std::vector<double> lambda_poly = multiply_polynomials(
		multiply_polynomials({ gram(0, 0), -1 }, 2, { gram(1, 1), -1 }, 2),
		3,
		{ gram(2, 2), -1 },
		2
	);

	lambda_poly[0] += gram(1, 0) * gram(2, 1) * gram(0 ,2) + gram(0, 1) * gram(1, 2) * gram(2, 0);

	std::vector<std::complex<double>> eigenvalues;
	Cardano_real_number(lambda_poly, eigenvalues);
	*/

	// Jacobi
	Matrix test(matrix_double{{ 1, 1, 3 }, { 1,5,1 }, { 3,1,1 }});
	//Matrix test(matrix_double{ { 4, -2, 2 }, { -2,5,0 }, { 2,0,6 } });
	double epsilon_elem = 0.0001;
	std::vector<double> eigenvalues_jacobi;
	Matrix T, U, S, V;

	//Jacobi(test, 3, epsilon_elem, 20, eigenvalues_jacobi, T);

	//SVD_matrix(test, cols, rows, eigenvalues_jacobi, T, U, S, V);

	Matrix RES;
	RES = V * S * U.transposition();

	return 0;
}