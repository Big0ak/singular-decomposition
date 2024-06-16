#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "../math/math_utils.h"
#include "../math/matrix.h"

class Jacobi
{
public:
	Jacobi(Matrix _a, int _n);
	~Jacobi();

	/*
	* input:
	*   epsilon_elem - maximum non-diagonal element
	*   max_rotation - max count rotation
	*/
	bool solve(double epsilon_elem, int max_rotation);
	std::vector<double> get_eigenvalues();
	Matrix get_eigenvectors();
	int get_count_rotation();

private:
	// find maximum element modem under diagonal
	bool find_max(Matrix &a, int n, double epsilon, std::pair<int, int>& max_index);

	Matrix a; // symmetric square matrix
	int n; // size matrix 'a'
	std::vector<double> eigenvalues; // eigen values
	Matrix T; // eigen vectors (calls)

	int count_rotation = 0;
};