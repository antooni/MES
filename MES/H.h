#pragma once

#include "kwadratury.h"
#include "calka.h"
#include "grid.h"
#include "matrix.h"
#include "jakobian.h"

struct H {
	Matrix* macierz;

	H(int nrElementu, Element4 element, Grid grid);
	~H() {
		delete macierz;
	}
};

struct Hbc {
	Matrix* macierz;

	Hbc(int nrElementu, Element4 element, Grid grid);
	~Hbc() {
		delete macierz;
	}

};

struct P {
	Matrix* macierz;

	P(int nrElementu, Element4 element, Grid grid);
	~P() {
		delete macierz;
	}
};

enum Walls {
	BOTTOM,
	RIGHT,
	TOP,
	LEFT,
};