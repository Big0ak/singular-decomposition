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
	Matrix(int _rows, int _cols);
	~Matrix();

	double& operator () (int i, int j) { return matrix[i][j]; }	// get elem
	void operator()(int i, int j, double data) { matrix[i][j] = data; } // set elem
	Matrix operator* (const Matrix&);
	Matrix operator+ (const Matrix&);
	Matrix operator- (const Matrix&);

	friend std::ostream& operator<< (std::ostream& fo, const Matrix& M);

	Matrix& transposition();
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
	rows = _matrix.size();
	cols = _matrix[0].size();
	matrix = _matrix;
}

Matrix::Matrix(int _rows, int _cols)
{
	matrix.resize(_rows);
	for (int i = 0; i < matrix.size(); i++)
		matrix[i].resize(_cols);
	cols = _cols;
	rows = _rows;
}

Matrix::~Matrix() { }

Matrix Matrix::operator* (const Matrix& M)
{
	Matrix res(rows, M.cols);
	if (cols == M.rows)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < M.cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += this->matrix[i][k] * M.matrix[k][j];
				res.matrix[i][j] = sum;
			}
	}
	return res;
}

Matrix Matrix::operator+ (const Matrix& M)
{
	Matrix res(*this);
	if (rows == M.rows && cols == M.cols)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				res.matrix[i][j] += M.matrix[i][j];
	}
	return res;
}

Matrix Matrix::operator- (const Matrix& M)
{
	Matrix res(*this);
	if (rows == M.rows && cols == M.cols)
	{
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				res.matrix[i][j] -= M.matrix[i][j];
	}
	return res;
}

std::ostream& operator<< (std::ostream& fo, const Matrix& M)
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

Matrix& Matrix::transposition()
{
	std::vector<std::vector<double>> m_new;
	m_new.resize(cols);
	for (int i = 0; i < m_new.size(); i++)
		m_new[i].resize(rows);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			m_new[j][i] = matrix[i][j];
		}
	}
	matrix = m_new;
	std::swap(cols, rows);
	return *this;
}
