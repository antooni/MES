#include "C.h"

#include <chrono>

constexpr auto SPECIFIC_HEAT = 700.0;
constexpr auto CONDUCTIVITY = 7800.0;


using namespace std::chrono;

C::C(int nrElementu, Element4 element, Grid grid)
{
	Matrix* C = new Matrix(element.N_FUNKCJI, element.N_FUNKCJI, 0.0);

	Matrix* N = new Matrix(element.N_FUNKCJI, element.indeksyPunktow.size(), 0.0);

	//dla kazdego punktu calkowania
	for (int k = 0; k < element.indeksyPunktow.size(); k++) {
		Jakobian* jakobian = new Jakobian(nrElementu, k, element, grid);

		double ksi = element.punktyCalkowania[element.indeksyPunktow[k]->xIndex]->x;
		double eta = element.punktyCalkowania[element.indeksyPunktow[k]->yIndex]->x;

		//dla kazdej funkcji ksztaltu
		for (int i = 0; i < element.N_FUNKCJI; i++) {
			double nvalue = element.N[i](ksi, eta);
			N->A[k][i] = nvalue;
		}
		//N->print();

		for (int i = 0; i < element.N_FUNKCJI; i++) {
			for (int j = 0; j < element.N_FUNKCJI; j++) {
				double n1 = N->A[k][i];
				double n2 = N->A[k][j];

				double a = n1 * n2;

				double waga = element.punktyCalkowania[element.indeksyPunktow[k]->xIndex]->A;
				waga *= element.punktyCalkowania[element.indeksyPunktow[k]->yIndex]->A;


				C->A[i][j] += SPECIFIC_HEAT * CONDUCTIVITY * a * jakobian->det * waga ;
			}
		}
		//C->print();
		delete jakobian;
	}

	delete N;
	macierz = C;
}
