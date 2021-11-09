#pragma once

#include "kwadratury.h"
#include "calka.h"
#include "grid.h"
#include "matrix.h"
#include "jakobian.h"

struct H {
	Matrix* macierz;

	H(Jakobian jakobian, Pochodne element, Grid grid);
};