#pragma once

struct Matrix {
	double** a;
	int n;
	int m;

	Matrix();
	Matrix(int _n, int _m, double** _a);
	Matrix(int _n, int _m, double value);

	void print();

	Matrix* getTransposed();

	Matrix add(Matrix* m1);
};