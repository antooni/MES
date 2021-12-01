#pragma once
#include "grid.h"
#include "jakobian.h"

class C {
public:
	Matrix* macierz;

	C(int nrElementu, Element4 element, Grid grid);
};