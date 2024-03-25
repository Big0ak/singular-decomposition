#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <complex>

#include "matrix.h"

#define M_PI	3.14159265358979323846

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

		answer[1].real(-(alpha + beta) / 2.0);
		answer[1].imag(std::sqrt(3) * ((alpha - beta) / 2));

		answer[2].real(-(alpha + beta) / 2.0);
		answer[2].imag(-std::sqrt(3) * ((alpha - beta) / 2));
		
		return;
	}
	else if (D < 0)
	{
		double phi = std::atan(-2 * std::sqrtf(D) / q);
		double r = 2 * std::sqrt(-p / 3);

		answer[0].real(r * std::cos(phi / 3.0));
		answer[0].imag(0.);

		answer[1].real(r * std::cos(phi + 2 * M_PI / 3.0));
		answer[1].imag(0.);

		answer[2].real(r * std::cos(phi + 4 * M_PI / 3.0));
		answer[2].imag(0.);
	}
	else {
		double alpha = Cube_root((-q / 2.0));
	
		answer[0].real(2.0 * alpha);
		answer[0].imag(0.);

		answer[1].real(-alpha);
		answer[1].imag(0.);

		answer[2].real(-alpha);
		answer[2].imag(0.);
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

	//Matrix a_t(m);
	//a_t.transposition();

	//Matrix a(m);
	//Matrix gram(a_t * a);

	std::vector<double> lambda_poly = multiply_polynomials(
		multiply_polynomials({ m[0][0], -1 }, 2, { m[1][1], -1 }, 2),
		3,
		{ m[2][2], -1 },
		2
	);

	lambda_poly[0] += m[1][0] * m[2][1] * m[0][2] + m[0][1] * m[1][2] * m[2][0];

	std::vector<std::complex<double>> answer;
	Cardano_real_number(lambda_poly, answer);

	return 0;
}