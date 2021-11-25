#pragma once

#include "kwadratury.h"
#include "calka.h"
#include "grid.h"
#include "matrix.h"
#include "jakobian.h"

struct H {
	Matrix* macierz;

	H(int nrElementu, Element4 element, Grid grid);
};

struct Hbc {
	Matrix* macierz;

	Hbc(int nrElementu, Element4 element, Grid grid);

};

struct P {
	Matrix* macierz;

	P(int nrElementu, Element4 element, Grid grid);
};

enum Walls {
	BOTTOM,
	RIGHT,
	TOP,
	LEFT,
};