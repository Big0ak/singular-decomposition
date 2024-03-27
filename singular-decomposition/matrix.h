#pragma once
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
	Matrix(int _cols, int _rows);
	~Matrix();

	double& operator()(int i, int j) { return matrix[i][j]; }	// получение элемента
	//void operator()(int i, int j, double data) { matrix[i][j] = data; }   // установка элемента
	Matrix operator* (const Matrix _matrix);
	friend std::ostream& operator << (std::ostream& fo, const Matrix& M);

	void transposition();
	void multiplication();
	void reverse();
	void determinant();
};

Matrix::Matrix(const Matrix& M)
{
	matrix.resize(M.rows);
	for (int i = 0; i < matrix.size(); i++)
		matrix[i].resize(M.cols);
	
	rows = M.rows;
	cols = M.cols;

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			matrix[i][j] = M.matrix[i][j];
}

Matrix::Matrix(std::vector<std::vector<double>> _matrix)
{
	cols = _matrix.size();
	rows = _matrix[0].size();
	matrix = _matrix;
}

Matrix::Matrix(int _cols, int _rows)
{
	matrix.resize(_cols);
	for (int i = 0; i < matrix.size(); i++)
		matrix[i].resize(_rows);
	cols = _cols;
	rows = _rows;
}

Matrix::~Matrix()
{

}

Matrix Matrix::operator* (const Matrix _matrix)
{
	Matrix res(rows, _matrix.cols);
	if (cols == _matrix.rows)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < _matrix.cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += this->matrix[i][k] * _matrix.matrix[k][j];
				res.matrix[i][j] = sum;
			}
	}
	return res;
}

std::ostream& operator << (std::ostream& fo, const Matrix& M)
{
	for (int i = 0; i < M.rows; i++)
	{
		fo << "  ";
		for (int j = 0; j < M.cols; j++)
			fo << M.matrix[i][j] << " \t";
		fo << std::endl;
	}
	return fo;
}

void Matrix::transposition()
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
	matrix = m_new;
}
