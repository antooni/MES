#include "matrix.h"
#include <iostream>


Matrix::Matrix()
{
	nKolumn = NULL;
	mRzedow = NULL;
	A = NULL;
}

Matrix::Matrix(int _nKolumn, int _mRzedow, double** _a)
{
	nKolumn = _nKolumn;
	mRzedow = _mRzedow;
	A = _a;
}

Matrix::Matrix(int _nKolumn, int _mRzedow, double value)
{
	nKolumn = _nKolumn;
	mRzedow = _mRzedow;

	A = new double* [mRzedow];
	for (int i = 0; i < mRzedow; i++) {
		A[i] = new double[nKolumn];
		for (int j = 0; j < nKolumn; j++) {
			A[i][j] = value;
		}
	}
}

void Matrix::print()
{
	std::cout << std::endl;
	for (int i = 0; i < mRzedow; i++) {
		for (int j = 0; j < nKolumn; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

Matrix* Matrix::getTransposed()
{
	Matrix* T = new Matrix(mRzedow,nKolumn,0.0);

	for (int i = 0; i < mRzedow; i++) {
		for (int j = 0; j < nKolumn; j++) {
			T->A[j][i] = A[i][j];
		}
	}

	return T;
}

Matrix Matrix::add(Matrix* m1)
{
	Matrix* result = new Matrix(mRzedow, nKolumn, 0.0);

	for (int i = 0; i < mRzedow; i++) {
		for (int j = 0; j < nKolumn; j++) {
			result->A[i][j] = A[i][j] + m1->A[i][j];
		}
	}

	return *result;
}





