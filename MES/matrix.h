#pragma once

struct Matrix {
	double** A;
	int nKolumn;
	int mRzedow;

	Matrix();
	~Matrix();
	Matrix(int _nKolumn, int _mRzedow, double** _A);
	Matrix(int _nKolumn, int _mRzedow, double value);

	void print();

	Matrix* getTransposed();

	Matrix* add(Matrix* m1);
	double min();
	double max();
};