#include "C.h"
#include <chrono>

constexpr auto SPECIFIC_HEAT = 700.0;
constexpr auto CONDUCTIVITY = 7800.0;

using namespace std::chrono;

C::C(int nrElementu, Element4 element, Grid grid)
{
	int m_punktow_calkowania = element.indeksyPunktow.size();

	Matrix* C = new Matrix(element.N_FUNKCJI, element.N_FUNKCJI, 0.0);

	Matrix* N = new Matrix(element.N_FUNKCJI, element.indeksyPunktow.size(), 0.0);

	//dla kazdego punktu calkowania
	for (int k = 0; k < m_punktow_calkowania; k++) {
		//wyliczamy jakobian dla kazdego punktu calkowania
		Jakobian* jakobian = new Jakobian(nrElementu, k, element, grid);

		//pobieramy wartosci ksi i eta w tym punkcie calkowania
		double ksi = element.punktyCalkowania[element.indeksyPunktow[k]->xIndex]->x;
		double eta = element.punktyCalkowania[element.indeksyPunktow[k]->yIndex]->x;

		//wyliczamy wartosci kolejnych funkcji kszta³tu (rz¹d)
		for (int i = 0; i < element.N_FUNKCJI; i++) {
			N->A[k][i] = element.N[i](ksi, eta);
		}

		//macierz przemnozona z transponowana i wspolczynnikami
		//dla kazdej komorki macierzy
		for (int i = 0; i < element.N_FUNKCJI; i++) {
			for (int j = 0; j < element.N_FUNKCJI; j++) {
				double n1 = N->A[k][i];
				double n2 = N->A[k][j];
				//transponowana komorka
				double T_komorka = n1 * n2;
				//waga
				double waga = element.punktyCalkowania[element.indeksyPunktow[k]->xIndex]->A;
				waga *= element.punktyCalkowania[element.indeksyPunktow[k]->yIndex]->A;

				//calkowanie
				C->A[i][j] += SPECIFIC_HEAT * CONDUCTIVITY * T_komorka * jakobian->det * waga ;
			}
		}
		delete jakobian;
	}

	delete N;
	macierz = C;
}
