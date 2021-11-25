#include "matrix.h"
#include <iostream>


Matrix::Matrix()
{
	n = NULL;
	m = NULL;
	a = NULL;
}

Matrix::Matrix(int _n, int _m, double** _a)
{
	n = _n;
	m = _m;
	a = _a;
}

Matrix::Matrix(int _n, int _m, double value)
{
	n = _n;
	m = _m;

	a = new double* [m];
	for (int i = 0; i < m; i++) {
		a[i] = new double[n];
		for (int j = 0; j < n; j++) {
			a[i][j] = value;
		}
	}
}

void Matrix::print()
{
	std::cout << std::endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << a[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

Matrix* Matrix::getTransposed()
{
	Matrix* T = new Matrix(m,n,0.0);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			T->a[j][i] = a[i][j];
		}
	}

	return T;
}

Matrix Matrix::add(Matrix* m1)
{
	Matrix* result = new Matrix(n, m, 0.0);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result->a[i][j] = a[i][j] + m1->a[i][j];
		}
	}

	return *result;
}





