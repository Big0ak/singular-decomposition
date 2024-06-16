#pragma once

#include <vector>
#include <complex>

#include "math_utils.h"

// writing polynomial as [0, x, x^2, x^3 x^4]
namespace Equation
{
	/*
	* input:
	*	poly - polynomial of the secont degree
	* result:
	*	answer - two(one) roots
	* return: count roots
	*/
	inline int Cubic_equations(std::vector<double> poly, std::vector<std::complex<double>>& answer)
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
	inline bool Cardano_real_number(std::vector<double> poly, std::vector<std::complex<double>>& answer)
	{
		if (poly.size() != 4)
		{
			return false;
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

		return true;
	}

	/*
	* input:
	*	poly - polynomial of the fourth degree
	* result:
	*	answer - four roots
	*/
	inline bool Ferrari_real_number(std::vector<double> poly, std::vector<std::complex<double>>& answer)
	{
		if (poly.size() != 5)
		{
			return false;
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

		return true;
	}
}
