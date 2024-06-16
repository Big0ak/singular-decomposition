#pragma once

#include <iostream>
#include <vector>

class Matrix
{
protected:
	std::vector<std::vector<double>> matrix;
	int cols, rows;
public:
	Matrix() : cols(0), rows(0) {}
	Matrix(const Matrix& M);
	Matrix(std::vector<std::vector<double>> _matrix);
	Matrix(int _rows, int _cols);
	~Matrix();

	double& operator () (int i, int j);	// get elem
	void operator()(int i, int j, double data); // set elem
	Matrix operator* (const Matrix&);
	Matrix operator+ (const Matrix&);
	Matrix operator- (const Matrix&);

	friend std::ostream& operator<< (std::ostream& fo, const Matrix& M);

	Matrix& transposition();
	void reverse();
	void determinant();
};