#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <complex>
#include <algorithm>

#include "matrix.h"

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

double Cube_root(double x)
{
	if (x > 0)
	{
		return std::pow(x, 1.0 / 3.0);
	}
	else {
		return -std::pow(-x, 1.0 / 3.0);
	}
}

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
	* x^3 + px + q = 0;
	*/

	double p = poly[1] - std::pow(poly[2], 2) / 3.0;
	double q = 2 * std::pow(poly[2], 3) / 27.0 - (poly[1] * poly[2]) / 3.0 + poly[0];

	double D = std::pow(q, 2) / 4.0 + std::pow(p, 3) / 27.0;
	if (D > 0)
	{
		// x0 - R, x1, x2 - C
		double alpha = Cube_root((-q / 2.0 + std::sqrt(D)));
		double beta = Cube_root((-q / 2.0 - std::sqrt(D)));

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
		double alpha = Cube_root((-q / 2.0));
	
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

// find maximum element model under diagonal
bool Find_max(Matrix& a, int n, double epsilon, std::pair<int, int>& max_index)
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

double sign(double x)
{
	if (x < 0)
		return -1.;
	
	if (x > 0)
		return 1.;

	return 0;
}

/*
* input:
*	m - symmetric square matrix
*   epsilon_elem - maximum non-diagonal element
*   max_rotation - max count rotation
* 
* result:
*	
*/
void Jacobi(Matrix a, int n, double epsilon_elem, int max_rotation, std::vector<double>& eigenvalues, Matrix& T)
{
	bool flag = true;
	int count_rotation = 0;
	double key_elem = 0;
	std::vector<double> b_im;
	b_im.resize(n);
	std::vector<double> b_jm;
	b_jm.resize(n);
	
	// work only with under diagonal elements
	while (flag)
	{
		std::pair<int, int> index;
		if (Find_max(a, n, epsilon_elem, index))
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
				s = std::sqrt(0.5 - r * sign(p * q));
			}

			// multiplication: B = T'*A*T
			
			// diagonal elements 
			double b_ii = std::pow(c, 2) * a(index.first, index.first) + std::pow(s, 2) * a(index.second, index.second) + 2 * c * s * a(index.first, index.second);
			double b_jj = std::pow(c, 2) * a(index.first, index.first) + std::pow(s, 2) * a(index.second, index.second) - 2 * c * s * a(index.first, index.second);
			a(index.first, index.first, b_ii);
			a(index.second, index.second, b_jj);

			a(index.first, index.second, 0); // i,j
			a(index.second, index.first, 0); // j,i

			// non-diagonal elements
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
				if (m != index.first && m != index.second)
				{
					a(index.first, m, b_im[m]);
					a(m, index.first, b_im[m]);

					a(index.second, m, b_jm[m]);
					a(m, index.second, b_jm[m]);
				}
			}

			// calculation eigenvectors
			if (count_rotation == 0)
			{
				T = Matrix(n, n);
				T(index.first, index.first, c);
				T(index.second, index.second, c);
				T(index.first, index.second, -s);
				T(index.second, index.first, s);
			}
			else {
				// T = T * T_new;
			}
		}
		else {
			flag = false;
		}
		
		count_rotation++;
		if (count_rotation >= max_rotation)
			flag = false;
	}

	eigenvalues.resize(n);
	for (int i = 0; i < n; i++)
		eigenvalues[i] = a(i, i);
}

void SVD_value(Matrix a, std::vector<double> eigenvalues, Matrix eigenvector, Matrix& U, Matrix& S, Matrix& V)
{
	std::sort(eigenvalues.begin(), eigenvalues.end());
	double r = -1;
	for (int i = eigenvalues.size() - 1; i >= 0; i--)
	{
		if (eigenvalues[i] == 0)
			r = i;
	}


}

int main()
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
	double epsilon_elem = 0.0001;
	std::vector<double> eigenvalues;
	Matrix T;

	Jacobi(test, 3, epsilon_elem, 20, eigenvalues, T);

	//SVD_value();

	return 0;
}