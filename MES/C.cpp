#include "C.h"

C::C(int nrElementu, Element4 element, Grid grid)
{
	Matrix* C = new Matrix(4, 4, 0.0);

	Matrix* N = new Matrix(4, 4, 0.0);

	for (int i = 0; i < 4; i++) {
		double ksi = element.kwadratura->stopnie[0]->punkty[element.indeksyPunktow[i]->xIndex]->x;
		double eta = element.kwadratura->stopnie[0]->punkty[element.indeksyPunktow[i]->yIndex]->x;
		for (int j = 0; j < 4; j++) {
			double nvalue = element.N[j](ksi, eta);
			N->A[i][j] = nvalue;
		}
	}

	double detJ = 0.00015625;

	for (int k = 0; k < 4; k++) {

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {

				Jakobian* jakobian = new Jakobian(nrElementu, j, element, grid);
				double n1 = N->A[k][i];
				double n2 = N->A[k][j];

				double a = n1 * n2;


				C->A[i][j] += 700.0 * 7800.0 * a * jakobian->det;
			}
		}
	}
	macierz = C;
}
