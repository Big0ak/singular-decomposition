#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <complex>
#include <algorithm>
#include <cstring>

#include "matrix.h"
#include "math_utils.h"

#define M_PI	3.14159265358979323846
#define matrix_double std::vector<std::vector<double>>

void print(std::vector<double> a, int n)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << a[i] << " ";
	}
	std::cout << std::endl;
}

// [0, x, x^2, x^3 ...]
std::vector<double> multiply_polynomials(std::vector<double> a, int n, std::vector<double> b, int m)
{
	int size = m + n - 1;
	std::vector<double> res;
	res.resize(size, 0);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			res[i + j] += a[i] * b[j];
			//print(res, size);
		}
	}

	return res;
}

/*
* input:
*	poly - polynomial of the secont degree
* result:
*	answer - two roots
* return: count roots
*/
int Cubic_equations(std::vector<double> poly, std::vector<std::complex<double>>& answer)
{
	if (poly.size() != 3)
	{
		return 0;
	}

	int count_root = 2;
	double D = std::pow(poly[1], 2) - 4 * poly[2] * poly[0];
	double D_sqrt = std::sqrt(std::abs(D));
	if (D > 0)
	{
		answer.resize(2);
		answer[0].real((-poly[1] + D_sqrt) / (2. * poly[2]));
		answer[0].imag(0.);

		answer[1].real((-poly[1] - D_sqrt) / (2. * poly[2]));
		answer[1].imag(0.);
	}
	else if (D < 0)
	{
		answer.resize(2);
		answer[0].real(-poly[1] / (2. * poly[2]));
		answer[0].imag(D_sqrt / (2. * poly[2]));

		answer[1].real(-poly[1] / (2. * poly[2]));
		answer[1].imag(-D_sqrt / (2. * poly[2]));
	}
	else {
		answer.resize(1);
		answer[0].real(-poly[1] / (2. * poly[2]));
		answer[0].imag(0.);
		count_root = 1;
	}

	return count_root;
}

/*
* input:
*	poly - polynomial of the third degree
* result:
*	answer - three roots
*/
void Cardano_real_number(std::vector<double> poly, std::vector<std::complex<double>>& answer)
{
	if (poly.size() != 4)
	{
		return;
	}
	answer.resize(3);

	for (int i = 2; i >= 0; i--)
	{
		poly[i] /= poly[3];
	}
	poly[3] = 1;

	/*
	* kardano
	* x^3 + px + q = 0
	*/

	double p = poly[1] - std::pow(poly[2], 2) / 3.0;
	double q = 2 * std::pow(poly[2], 3) / 27.0 - (poly[1] * poly[2]) / 3.0 + poly[0];

	double D = std::pow(q, 2) / 4.0 + std::pow(p, 3) / 27.0;
	if (D > 0)
	{
		// x0 - R, x1, x2 - C
		double alpha = MATH_utils::cube_root((-q / 2.0 + std::sqrt(D)));
		double beta = MATH_utils::cube_root((-q / 2.0 - std::sqrt(D)));

		answer[0].real(alpha + beta);
		answer[0].imag(0.);

		answer[1].real(-(alpha + beta) / 2.0);
		answer[1].imag(std::sqrt(3) * ((alpha - beta) / 2));

		answer[2].real(-(alpha + beta) / 2.0);
		answer[2].imag(-std::sqrt(3) * ((alpha - beta) / 2));
	}
	else if (D < 0)
	{
		// x0, x1, x2 - R
		std::complex<double> t;
		t.real(-q / 2);
		t.imag(std::sqrt(-D));
		double t_r = std::abs(t);
		double t_phi = std::arg(t);

		double phi = 0;
		if (q < 0)
			phi = std::atan(-2 * std::sqrt(-D) / q);
		else if (q > 0)
			phi = std::atan(-2 * std::sqrt(-D) / q) + M_PI;
		else
			phi = M_PI / 2;
			
		//double r = std::sqrt(-std::pow(p, 3) / 27);
		double r = 2 * std::sqrt(-p / 3);

		answer[0].real(r * std::cos(phi / 3.0));
		answer[0].imag(0.);

		answer[1].real(r * std::cos((phi + 2 * M_PI) / 3.0));
		answer[1].imag(0.);

		answer[2].real(r * std::cos((phi + 4 * M_PI) / 3.0));
		answer[2].imag(0.);
	}
	else {
		// x0, x1 == x2 - R
		double alpha = MATH_utils::cube_root((-q / 2.0));
	
		answer[0].real(2.0 * alpha);
		answer[0].imag(0.);

		answer[1].real(-alpha);
		answer[1].imag(0.);

		answer[2].real(-alpha);
		answer[2].imag(0.);
	}

	answer[0].real(answer[0].real() - poly[2] / (3 * poly[3]));
	answer[1].real(answer[1].real() - poly[2] / (3 * poly[3]));
	answer[2].real(answer[2].real() - poly[2] / (3 * poly[3]));

	return;
}

/*
* input:
*	poly - polynomial of the fourth degree
* result:
*	answer - four roots
*/
void Ferrari_real_number(std::vector<double> poly, std::vector<std::complex<double>>& answer)
{
	if (poly.size() != 5)
	{
		return;
	}
	answer.resize(4);

	for (int i = 3; i >= 0; i--)
	{
		poly[i] /= poly[4];
	}
	poly[4] = 1;

	/*
	* calculate cubic resolvent
	* 2s^3 - ps^2 - 2rs + rp - q^4/4 = 0
	*/ 
	std::vector<double> resolvent;
	std::vector<std::complex<double>> answer_resolvent;
	resolvent.resize(4);

	double p, q, r;
	p = poly[2] - 3 * std::pow(poly[3], 2) / 8;
	q = std::pow(poly[3], 3) / 8 - poly[3] * poly[2] / 2 + poly[1];
	r = poly[0] - poly[1] * poly[3] / 4 + std::pow(poly[3], 2) * poly[2] / 16 - 3 * std::pow(poly[3], 4) / 256;
	resolvent[0] = r * p - std::pow(q, 2) / 4;
	resolvent[1] = -2 * r;
	resolvent[2] = -p;
	resolvent[3] = 2;

	Cardano_real_number(resolvent, answer_resolvent);

	double real_ans = 0;
	for (int i = 0; i < 3; i++)
	{
		if (answer_resolvent[i].imag() == 0)
		{
			real_ans = answer_resolvent[i].real();
			break;
		}
	}

	// solving cubic equations
	std::vector<std::complex<double>> answer_cubic;
	std::vector<double> poly_cubic;
	poly_cubic.resize(3);
	poly_cubic[0] = q / (2 * std::sqrt(2 * real_ans - p)) + real_ans;
	poly_cubic[1] = -std::sqrt(2 * real_ans - p);
	poly_cubic[2] = 1.;
	if (Cubic_equations(poly_cubic, answer_cubic) == 2) // !!! == 2 ???
	{
		answer[0] = answer_cubic[0];
		answer[1] = answer_cubic[1];
	}

	answer_cubic.clear();
	poly_cubic[0] = -q / (2 * std::sqrt(2 * real_ans - p)) + real_ans;
	poly_cubic[1] = std::sqrt(2 * real_ans - p);
	poly_cubic[2] = 1.;
	if (Cubic_equations(poly_cubic, answer_cubic) == 2)
	{
		answer[2] = answer_cubic[0];
		answer[3] = answer_cubic[1];
	}

	answer[0].real(answer[0].real() - poly[3] / (4 * poly[4]));
	answer[1].real(answer[1].real() - poly[3] / (4 * poly[4]));
	answer[2].real(answer[2].real() - poly[3] / (4 * poly[4]));
	answer[3].real(answer[3].real() - poly[3] / (4 * poly[4]));
}

// find maximum element modem under diagonal
bool find_max(Matrix& a, int n, double epsilon, std::pair<int, int>& max_index)
{
	max_index.first = -1;
	max_index.second = -1;
	double max = -1;
	for (int i = 1; i < n; i++)
		for (int j = 0; j < i; j++)
			if (std::abs(a(i, j)) > max && std::abs(a(i, j)) > epsilon)
			{
				max = std::abs(a(i, j));
				max_index.first = i;
				max_index.second = j;
			}
	return max == -1 ? false : true;
}

/*
* input:
*	m - symmetric square matrix
*   epsilon_elem - maximum non-diagonal element
*   max_rotation - max count rotation
* 
* result:
*	eigenvalues - eigen values
*	T - eigen vectors (calls)
*/
void Jacobi(Matrix a, int n, double epsilon_elem, int max_rotation, std::vector<double>& eigenvalues, Matrix& T)
{
	bool flag = true;
	int count_rotation = 0;
	double key_elem = 0;
	std::vector<double> b_im, b_jm, t_mi, t_mj;
	b_im.resize(n);
	b_jm.resize(n);
	t_mi.resize(n);
	t_mj.resize(n);
	
	// work only with under diagonal elements
	while (flag)
	{
		std::pair<int, int> index;
		if (find_max(a, n, epsilon_elem, index))
		{
			key_elem = a(index.first, index.second);
			double p = 2 * key_elem;
			double q = a(index.first, index.first) - a(index.second, index.second);
			double d = std::sqrt(std::pow(p, 2) + std::powf(q, 2));

			/*
			* rotation matrix
			* T = |c -s|
			*	  |s  c|
			*/ 			
			double c, s;
			if (q == 0)
			{
				c = std::sqrt(2) / 2.;
				s = c;
			}
			else {
				double r = std::abs(q) / (2. * d);
				c = std::sqrt(0.5 + r);
				s = std::sqrt(0.5 - r) * MATH_utils::sign(p * q);
			}

			// multiplication: B = T'*A*T
			std::fill(b_im.begin(), b_im.end(), 0);
			std::fill(b_jm.begin(), b_jm.end(), 0);

			// (i,j) = (j,i) = 0

			// 1. diagonal elements 
			b_im[index.first] = std::pow(c, 2) * a(index.first, index.first) + std::pow(s, 2) * a(index.second, index.second) + 2 * c * s * a(index.first, index.second);
			b_jm[index.second] = std::pow(c, 2) * a(index.first, index.first) + std::pow(s, 2) * a(index.second, index.second) - 2 * c * s * a(index.first, index.second);

			// 2. non-diagonal elements
			// calculate elements
			for (int m = 0; m < n; m++)
			{
				if (m != index.first && m != index.second)
				{
					b_im[m] = c * a(m, index.first) + s * a(m, index.second);
					b_jm[m] = -s * a(m, index.first) + c * a(m, index.second);
				}
			}

			// set elements
			for (int m = 0; m < n; m++)
			{
				a(index.first, m, b_im[m]);
				a(m, index.first, b_im[m]);

				a(index.second, m, b_jm[m]);
				a(m, index.second, b_jm[m]);
			}

			// calculation eigenvectors
			if (count_rotation == 0)
			{
				T = Matrix(n, n);
				for (int i = 0; i < n; i++)
					T(i, i, 1);

				T(index.first, index.first, c);
				T(index.second, index.second, c);
				T(index.first, index.second, -s);
				T(index.second, index.first, s);
			}
			else {
				// T = T * T_new;
				std::fill(t_mi.begin(), t_mi.end(), 0);
				std::fill(t_mj.begin(), t_mj.end(), 0);
				// calculate elements
				for (int m = 0; m < n; m++)
				{
					t_mi[m] = c * T(m, index.first) + s * T(m, index.second);
					t_mj[m] = -s * T(m, index.first) + c * T(m, index.second);
				}

				// set elements
				for (int m = 0; m < n; m++)
				{
					T(m, index.first, t_mi[m]);
					T(m, index.second, t_mj[m]);
				}
			}

			count_rotation++;
			if (count_rotation >= max_rotation)
				flag = false;
		}
		else {
			flag = false;
		}
	}

	eigenvalues.resize(n);
	for (int i = 0; i < n; i++)
		eigenvalues[i] = a(i, i);
}

void normalization(Matrix& eigenvectors, int n)
{
	for (int j = 0; j < n; j++)
	{
		double norm = 0;
		for (int i = 0; i < n; i++)
		{
			double tmp = eigenvectors(i, j);
			norm += tmp * tmp;
		}
		norm = std::sqrt(norm);
		for (int i = 0; i < n; i++)
		{
			eigenvectors(i, j, eigenvectors(i, j) / norm);
		}
	}
}

void SVD_matrix(Matrix A, int n, int m, std::vector<double> eigenvalues, Matrix eigenvectors, Matrix& U, Matrix& S, Matrix& V)
{
	// n > m
	U = eigenvectors; // m * m
	normalization(U, m);

	std::vector<double> sindular_num = eigenvalues; // 1 * m
	for (int i = 0; i < m; i++)
		sindular_num[i] = std::sqrt(sindular_num[i]); // negative eigen values not possible
	std::sort(sindular_num.begin(), sindular_num.end(), std::greater<>());

	// find index zero singular number
	int r = m;
	int i = r - 1;
	while (i < m && i == r - 1)
	{
		if (sindular_num[i] == 0)
			r = i;
		i--;
	}

	// calculation V - n * n
	V = A * U; // n * m
	// add to n*n
	for (int j = 0; j < r; j++)
	{
		double sig = sindular_num[j];
		for (int i = 0; i < n; i++)
		{
			V(i, j, V(i, j) / sig);
		}
	}

	// V add to ortogonal base (r+1...n) 

	S = Matrix(n, m); // n * m
	for (int i = 0; i < r; i++)
	{
		S(i, i, sindular_num[i]);
	}
}

void SVD(Matrix A, int n, int m, Matrix& U, Matrix& S, Matrix& V)
{
	//bool conjugate = false;
	if (m > n)
	{
		A.transposition();
		std::swap(n, m);
		//conjugate = true;
	}

	Matrix A_t(A);
	std::vector<double> eigenvalues;
	Matrix eigenvectors;

	Matrix gram;
	gram = A_t.transposition() * A; // n * n

	switch (n)
	{
		//case 3:
		//{
		//	//Cardano_real_number()
		//  // find eigenvectors
		//	break;
		//}
		//case 4:
		//{

		//	break;
		//}
		default:
		{
			double epsilon_elem = 0.0001;
			int max_iteration = 20;
			Jacobi(gram, m, epsilon_elem, 20, eigenvalues, eigenvectors);
			break;
		}
	}

	SVD_matrix(A, n, m, eigenvalues, eigenvectors, U, S, V);

	//if (conjugate)
	//{
	//	U.transposition();
	//	S.transposition();
	//	V.transposition();
	//}
}

// test svd
int main2()
{
	int n = 0, m = 0;

	//n = 2; m = 4;
	//Matrix A ({ { 1, 1, 1, 1 }, {1, 1, 1, 1} });

	n = 3; m = 3;
	Matrix A({ {1,1,3}, {1,5,1}, {3,1,1} });
	//Matrix A({ { 4, 3, 2 }, { -2,5,0 }, { 2,0,6 } });

	Matrix T, U, S, V;
	SVD(A, n, m, U, S, V);

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


//void random_matrix(double**& matrix, double*& xlist, double*& eigenvalues, int& minLamInd, int N, int range)
//{
//	eigenvalues = new double[N]; //соб значения матрицы
//	double* w = new double[N];
//	double* e = new double[N]; // единичный вектор
//	matrix = new double* [N];
//	for (int i = 0; i < N; i++)
//		matrix[i] = new double[i];
//	double	wlength = 0;
//
//	srand(time(0));
//	for (int i = 0; i < N; i++)
//	{
//		eigenvalues[i] = -range + 2 * (rand() % range);
//		if (fabs(eigenvalues[i]) < fabs(eigenvalues[minLamInd]))
//			minLamInd = i;
//		w[i] = -range + 2 * (rand() % range);
//		e[i] = 1;
//	}
//
//	wlength = 1 / norm(w, N);
//	for (int i = 0; i < N; i++)
//		w[i] *= wlength;
//	Matrix Lambda(N, N, eigenvalues, true); //диагональная матрица с соб значениями на диагнонали
//	Matrix W(N, 1, w); //вектор w
//	//cout << W;
//	Matrix Wt(1, N, w); //транспонированный вектор
//	Matrix E(N, N, e, true); //единичная матрица
//	Matrix H = E - (W * Wt) * 2;
//	//cout << H;
//	Matrix A = H * Lambda * H;
//	std::cout << A << std::endl;
//
//	xlist = new double[N];
//	for (int i = 0; i < N; i++)
//		xlist[i] = H(i, minLamInd); //компоненты реального собственного вектора, отвечающего макс. соб значению
//
//	for (int i = 0; i < N; i++)
//		for (int j = 0; j < N; j++)
//			matrix[i][j] = A(i, j);
//}


// Ferrari test
int main()
{
	std::vector<std::complex<double>> answer;
	std::vector<double> poly;
	poly.resize(5);
	poly[0] = -10;
	poly[1] = 16;
	poly[2] = -5;
	poly[3] = 3;
	poly[4] = 2;

	Ferrari_real_number(poly, answer);

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

	Jacobi(test, 3, epsilon_elem, 20, eigenvalues_jacobi, T);

	SVD_matrix(test, cols, rows, eigenvalues_jacobi, T, U, S, V);

	Matrix RES;
	RES = V * S * U.transposition();

	return 0;
}