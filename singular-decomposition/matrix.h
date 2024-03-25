#pragma once
#include <vector>

class Matrix
{
public:
	Matrix() : cols(0), rows(0) {}
	Matrix(std::vector<std::vector<double>> _matrix)
	{
		cols = _matrix.size();
		rows = _matrix[0].size();
		matrix = _matrix;
	}
	Matrix(int _cols, int _rows)
	{
		matrix.resize(_cols);
		for (int i = 0; i < matrix.size(); i++)
			matrix[i].resize(_rows);
		cols = _cols;
		rows = _rows;
	}


	~Matrix();

	double& operator()(int i, int j) { return matrix[i][j]; }	// получение элемента
	void operator()(int i, int j, double data) { matrix[i][j] = data; }   // установка элемента

	void transposition()
	{
		std::vector<std::vector<double>> m_new;
		m_new.resize(rows);
		for (int i = 0; i < m_new.size(); i++)
			m_new[i].resize(cols);
		for (int i = 0; i < cols; i++)
		{
			for (int j = 0; j < rows; j++)
			{
				m_new[j][i] = matrix[i][j];
			}
		}
		m_new = matrix;
	}

	Matrix& operator* (const Matrix _matrix)
	{
		
	}

	void multiplication();
	void reverse();
	void determinant();

private:
	std::vector<std::vector<double>> matrix;
	int cols, rows;
};