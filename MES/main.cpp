#include <iostream>
#include <chrono>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"
#include "jakobian.h"
#include "matrix.h"
#include "H.h"
#include "C.h"

//paraview, pliki txt

using namespace std::chrono;
const double eps = 1e-12;

//const double DELTA_T = 50.0;
const double DELTA_T = 1.0;

Matrix* gauss(Matrix* A, Matrix* B);

int main()
{
	Element4* uniwersalnyElement = new Element4(TRZY_PUNKTOWY);

	//Grid* g = new Grid(0.1, 0.1, 4, 4);
	Grid* g = new Grid(0.1, 0.1, 31, 31);
	//Grid* g = new Grid(0.25, 0.25, 2, 2);

	//g->homework();

	const int ILOSC_WEZLOW = g->nodes.size();
	
	//for (int iteracja = 0; iteracja < 10; iteracja++) {
	for (int iteracja = 0; iteracja < 3; iteracja++) {
		Matrix* HEX = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);
		Matrix* CEX = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);
		Matrix* PEX = new Matrix(1, ILOSC_WEZLOW, 0.0);

		for (int nrElementu = 0; nrElementu < g->elements.size(); nrElementu++) {
			H* h = new H(nrElementu, *uniwersalnyElement, *g);
			Hbc* hbc = new Hbc(nrElementu, *uniwersalnyElement, *g);
			P* p = new P(nrElementu, *uniwersalnyElement, *g);
			C* c = new C(nrElementu, *uniwersalnyElement, *g);

			//h->macierz->print();
			//hbc->macierz->print();
			//p->macierz->print();
			//c->macierz->print();

			//c->macierz->print();

			Matrix* res = h->macierz->add(hbc->macierz);

			for (int i = 0; i < res->mRzedow; i++) {
				for (int j = 0; j < res->nKolumn; j++) {
					HEX->A[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += res->A[i][j];
				}
			}

			for (int i = 0; i < res->mRzedow; i++) {
				for (int j = 0; j < res->nKolumn; j++) {
					CEX->A[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += c->macierz->A[i][j];
				}
			}

			for (int j = 0; j < p->macierz->mRzedow; j++) {
				PEX->A[g->elements[nrElementu]->nodesID[j] - 1][0] += p->macierz->A[j][0];
			}

			delete h;
			delete hbc;
			delete p;
			delete c;
		}
		//HEX->print();
		//CEX->print();
		//PEX->print();


		//lewa strona rownania
		Matrix* RES = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);

		for (int i = 0; i < ILOSC_WEZLOW; i++)
			for (int j = 0; j < ILOSC_WEZLOW; j++) {

				RES->A[i][j] = HEX->A[i][j] + (CEX->A[i][j] / DELTA_T);
			}

		//prawa strona rownania
		Matrix* PES = new Matrix(1, ILOSC_WEZLOW, 0.0);

		for (int i = 0; i < ILOSC_WEZLOW; i++)
			for (int j = 0; j < ILOSC_WEZLOW; j++) {
				PES->A[i][0] += g->nodes[j]->t0 * (CEX->A[i][j] / DELTA_T);
			}

		for (int i = 0; i < ILOSC_WEZLOW; i++) {
			PES->A[i][0] += PEX->A[i][0];
		}

		Matrix* X = gauss(RES, PES);
		std::cout << (iteracja + 1 ) * DELTA_T << " : " << X->min() << " " << X->max()<< std::endl;

		for (int i = 0; i < ILOSC_WEZLOW; i++) {
			double a = X->A[i][0];
			g->nodes[i]->t0 = X->A[i][0];
		}
		delete HEX;
		delete CEX;
		delete PEX;
		delete RES;
		delete PES;
		delete X;
	}
	return 0;
}

Matrix* gauss( Matrix* A, Matrix* B)
{
	int n = A->nKolumn;

	Matrix* AB = new Matrix(A->nKolumn+1, A->mRzedow, 0.0);

	for (int i = 0; i < AB->nKolumn; i++) {
		for (int j = 0; j < AB->mRzedow; j++) {
			if (i == AB->nKolumn-1) AB->A[j][i] = B->A[j][0];
			else AB->A[i][j] = A->A[i][j];
		}
	}

	int i, j, k;
	double m, s;

	Matrix* X = new Matrix(1, A->nKolumn, 0.0);

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB->A[i][i]) < eps) {
				std::cout << "DZIELNIK ZERO W GAUSS" << std::endl;
				return new Matrix(1, 1, 0.0);
			}
			m = -AB->A[j][i] / AB->A[i][i];
			for (k = i + 1; k <= n; k++)
				AB->A[j][k] += m * AB->A[i][k];
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		s = AB->A[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB->A[i][j] * X->A[j][0];
		if (fabs(AB->A[i][i]) < eps) {
			std::cout << "DZIELNIK ZERO W GAUSS" << std::endl;
			return new Matrix(1, 1, 0.0);
		}
		X->A[i][0] = s / AB->A[i][i];
	}

	delete AB;

	return X;
}