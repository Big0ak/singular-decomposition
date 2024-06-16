#include "Jacobi.h"

Jacobi::Jacobi(Matrix _a, int _n)
{
	a = _a;
	n = _n;
}

Jacobi::~Jacobi()
{

}

bool Jacobi::solve(double epsilon_elem, int max_rotation)
{
	count_rotation = 0;
	bool flag = true;
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
			double d = std::sqrt(std::pow(p, 2) + std::pow(q, 2));

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
			b_jm[index.second] = std::pow(s, 2) * a(index.first, index.first) + std::pow(c, 2) * a(index.second, index.second) - 2 * c * s * a(index.first, index.second);

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
			if (count_rotation > max_rotation)
				flag = false;
		}
		else {
			flag = false;
		}
	}

	if (count_rotation == max_rotation)
	{
		eigenvalues.clear();
		return false;
	}

	eigenvalues.resize(n);
	for (int i = 0; i < n; i++)
		eigenvalues[i] = a(i, i);

	return true;
}

// find maximum element modem under diagonal
bool Jacobi::find_max(Matrix& a, int n, double epsilon, std::pair<int, int>& max_index)
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

std::vector<double> Jacobi::get_eigenvalues()
{
	return eigenvalues;
}

Matrix Jacobi::get_eigenvectors()
{
	return T;
}

int Jacobi::get_count_rotation()
{
	return count_rotation;
}