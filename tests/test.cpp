#include "test.h"

void TESTS::jacobi_randomHouseholder_inaccuracy()
{

}

void TESTS::jacobi_severalRandomHouseholder_meanInaccuracy()
{

}

double norm(double* x, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		double tmp = x[i];
		sum += tmp * tmp;
	}
	return sqrt(sum);
}

//void random_matrix(double**& matrix, double*& xlist, double*& eigenvalues, int& minLamInd, int N, int range)
//{
//	eigenvalues = new double[N]; //��� �������� �������
//	double* w = new double[N];
//	double* e = new double[N]; // ��������� ������
//	matrix = new double* [N];
//	for (int i = 0; i < N; i++)
//		matrix[i] = new double[i];
//	double wlength = 0;
//
//	std::srand(time(0));
//	for (int i = 0; i < N; i++)
//	{
//		eigenvalues[i] = -range + 2 * (rand() % range);
//		if (fabs(eigenvalues[i]) < fabs(eigenvalues[minLamInd]))
//			minLamInd = i;
//		w[i] = -range + 2 * (rand() % range);
//		e[i] = 1;
//	}
//
//	wlength = 1 / norm(w, N);
//	for (int i = 0; i < N; i++)
//		w[i] *= wlength;
//	Matrix Lambda(N, N, eigenvalues, true); //������������ ������� � ��� ���������� �� ����������
//	Matrix W(N, 1, w); //������ w
//	//cout << W;
//	Matrix Wt(1, N, w); //����������������� ������
//	Matrix E(N, N, e, true); //��������� �������
//	Matrix H = E - (W * Wt) * 2;
//	//cout << H;
//	Matrix A = H * Lambda * H;
//	std::cout << A << std::endl;
//
//	xlist = new double[N];
//	for (int i = 0; i < N; i++)
//		xlist[i] = H(i, minLamInd); //���������� ��������� ������������ �������, ����������� ����. ��� ��������
//
//	for (int i = 0; i < N; i++)
//		for (int j = 0; j < N; j++)
//			matrix[i][j] = A(i, j);
//}