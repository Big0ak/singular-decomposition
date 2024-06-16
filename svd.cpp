#include "svd.h"

SVD::SVD()
{

}

SVD::SVD(Solving solving)
{

}

SVD::~SVD()
{

}

bool SVD::solve(Matrix A, int n, int m, Matrix& U, Matrix& S, Matrix& V)
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
		Jacobi jacobi(gram, m);
		if (jacobi.solve(epsilon_elem, 20) == false)
		{
			return false;
		}

		eigenvalues = jacobi.get_eigenvalues();
		eigenvectors = jacobi.get_eigenvectors();

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

void SVD::SVD_matrix(Matrix A, int n, int m, std::vector<double> eigenvalues, Matrix eigenvectors, Matrix& U, Matrix& S, Matrix& V)
{
	// n > m
	U = eigenvectors; // m * m
	normalization(U, m);

	std::vector<double> singular_num = eigenvalues; // 1 * m
	for (int i = 0; i < m; i++)
		singular_num[i] = std::sqrt(singular_num[i]); // negative eigen values not possible
	std::sort(singular_num.begin(), singular_num.end(), std::greater<>());

	// find index zero singular number
	int r = m;
	int i = r - 1;
	while (i < m && i == r - 1)
	{
		if (singular_num[i] == 0)
			r = i;
		i--;
	}

	// calculation V - n * n
	V = A * U; // n * m
	// add to n*n
	for (int j = 0; j < r; j++)
	{
		double sig = singular_num[j];
		for (int i = 0; i < n; i++)
		{
			V(i, j, V(i, j) / sig);
		}
	}

	// V add to ortogonal base (r+1...n) 

	S = Matrix(n, m); // n * m
	for (int i = 0; i < r; i++)
	{
		S(i, i, singular_num[i]);
	}
}

void SVD::normalization(Matrix& eigenvectors, int n)
{
	for (int j = 0; j < n; j++)
	{
		double norm = 0;
		for (int i = 0; i < n; i++)
		{
			double tmp = eigenvectors(i, j);
			norm += tmp * tmp;
		}
		norm = 1 / std::sqrt(norm);
		for (int i = 0; i < n; i++)
		{
			eigenvectors(i, j, eigenvectors(i, j) * norm);
		}
	}
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

//void print(std::vector<double> a, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		std::cout << a[i] << " ";
//	}
//	std::cout << std::endl;
//}