#pragma once
#include <vector>


class Matrix
{
public:
	Matrix() : cols(0), rows(0) {}

	~Matrix();
	void multiplication();
	void reverse();
	void transposition();
	void determinant();

private:
	std::vector<std::vector<double>> matrix;
	int cols, rows;
};