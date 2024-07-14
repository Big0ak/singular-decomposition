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
	bool conjugate = false;
	
	if (m > n)
	{
		A.transposition();
		std::swap(n, m);
		conjugate = true;
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
		int max_rotation = 50;
		Jacobi jacobi(gram, m);
		if (jacobi.solve(epsilon_elem, max_rotation) == false)
		{
			return false;
		}

		eigenvalues = jacobi.get_eigenvalues();
		eigenvectors = jacobi.get_eigenvectors();

		break;
	}
	}

	SVD_matrix(A, n, m, eigenvalues, eigenvectors, U, S, V);

	if (conjugate)
	{
		//U.transposition();
		//S.transposition(); // !!
		//V.transposition();
	}

	return true;
}

void SVD::SVD_matrix(Matrix A, int n, int m, std::vector<double> eigenvalues, Matrix eigenvectors, Matrix& U, Matrix& S, Matrix& V)
{
	// n - row, m - cols
	// necessary n > m 
	U = eigenvectors; // U - m * m
	normalization(U, m);

	std::vector<double> singular_num = eigenvalues; // 1 * m
	for (int i = 0; i < m; i++)
	{
		singular_num[i] = std::sqrt(singular_num[i]); // negative eigen values not possible
	}
		
	// sort singular num on decreasing
	for (int left = 0; left < m; left++)
		for (int i = m-1; i > left; i--)
		{
			if (singular_num[i] > singular_num[i - 1])
			{
				std::swap(singular_num[i], singular_num[i - 1]);
				U.swap_colums(i, i - 1);
			}
		}

	// find index zero singular number
	int r = m;
	int index = r - 1;
	while (index < m && index == r - 1)
	{
		if (singular_num[index] == 0)
			r = index;
		index--;
	}

	V = Matrix(n, n); // V - n * n
	for (int j = 0; j < r; j++) // r < m
	{
		// V[i] = A*U[j] / singular_num[i]
		double sig = singular_num[j];
		double x = 0;
		for (int i = 0; i < n; i++)
		{
			x = 0;
			for (int k = 0; k < m; k++)
			{
				x += A(i, k) * U(k, j);
			}
			V(i, j, x / sig);
		}
	}
	// have: V - n * r (r <= n <= m)
	// V add to ortogonal base (r+1...n)
	Addition_base(V, n, n, r);

	S = Matrix(n, m); // S - n * m
	for (int i = 0; i < r; i++)
	{
		S(i, i, singular_num[i]);
	}
}

bool SVD::Addition_base(Matrix& A, int n, int m, int r)
{
	if (r == m)
		return true;



	Matrix A_t(A);
	Matrix check = A_t * A; // if one matrix -> true;

	return true;
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