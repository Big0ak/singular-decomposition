#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include <ccomplex>

void print(std::vector<double> a, int n)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << a[i] << " ";
	}
	std::cout << std::endl;
}

// [0, x, x^2, x^3 ...]
std::vector<double> multiply_polynomials(std::vector<double> a, std::vector<double> b, int n, int m)
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

	double p = poly[1] - std::pow(poly[2], 2) / 3.0;
	double q = 2 * std::pow(poly[2], 3) / 27.0 - (poly[1] * poly[2]) / 3.0 + poly[0];

	double D = std::pow(q, 2) / 4.0 + std::pow(p, 3) / 27.0;
	if (D > 0)
	{
		// x0 - R, x1, x2 - C
		double alpha = std::pow((-q / 2.0 + std::sqrt(D)), 1.0 / 3.0);
		double beta = std::pow((-q / 2.0 - std::sqrt(D)), 1.0 / 3.0);

		answer[0].real(alpha + beta);
		answer[1].real(-(alpha + beta) / 2.0);
		answer[1].imag(std::sqrt(3) * ((alpha - beta) / 2));

		answer[2].real(-(alpha + beta) / 2.0);
		answer[2].imag(-std::sqrt(3) * ((alpha - beta) / 2));
		
		return;
	}
	else if (D < 0)
	{

	}
	else {
		double alpha = std::pow((-q / 2.0 + std::sqrt(D)), 1.0 / 3.0);
		double beta = std::pow((-q / 2.0 - std::sqrt(D)), 1.0 / 3.0);

		answer[0].real(2.0 * alpha);
		answer[0].imag(0);

		answer[1].real(-alpha);
		answer[1].imag(0);

		answer[2].real(-alpha);
		answer[2].imag(0);
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

	std::vector<double> lambda_poly = multiply_polynomials(
		multiply_polynomials({ m[0][0], -1 }, { m[1][1], -1 }, 2, 2),
		{ m[2][2], -1 },
		3,
		2
	);

	lambda_poly[0] += m[1][0] * m[2][1] * m[0][2] + m[0][1] * m[1][2] * m[2][0];

	/*
	 * kardano
	 * x^3 + px + q = 0;
	 */

	std::vector<std::complex<double>> answer;
	Cardano_real_number(lambda_poly, answer);


	return 0;
}