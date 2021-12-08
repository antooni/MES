#include <iostream>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"
#include "jakobian.h"
#include "matrix.h"
#include "H.h"
#include "C.h"

const double eps = 1e-12;

const double DELTA_T = 1.0;

double fun(double x) {return 5 * x * x + 3 * x + 6;}

double fun2(double x, double y) {return 5 * x * x * y * y + 3 * x * y + 6;}

Funkcja* funkcja = new Funkcja(fun, "f(x) = 5x^2 + 3x + 6");
Funkcja* funkcja2 = new Funkcja(fun2, "f(x) = 5x^2y^2 + 3xy + 6");

Matrix* gauss(Matrix* A, Matrix* B);


int main()
{
	Element4* element = new Element4(DWU_PUNKTOWY);

	Grid* g = new Grid(0.1, 0.1, 31, 31);

	//g->homework();

	const int ILOSC_WEZLOW = g->nodes.size();
	


	for (double iteracja = 0.0; iteracja < 2.0; iteracja += 1.0) {
		Matrix* HEX = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);
		Matrix* CEX = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);
		Matrix* PEX = new Matrix(1, ILOSC_WEZLOW, 0.0);

		for (int nrElementu = 0; nrElementu < g->elements.size(); nrElementu++) {
			H* h = new H(nrElementu, *element, *g);
			Hbc* hbc = new Hbc(nrElementu, *element, *g);
			P* p = new P(nrElementu, *element, *g);
			C* c = new C(nrElementu, *element, *g);

			Matrix res = h->macierz->add(hbc->macierz);


			for (int i = 0; i < res.mRzedow; i++) {
				for (int j = 0; j < res.nKolumn; j++) {
					HEX->A[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += res.A[i][j];
				}
			}

			for (int i = 0; i < res.mRzedow; i++) {
				for (int j = 0; j < res.nKolumn; j++) {
					CEX->A[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += c->macierz->A[i][j];
				}
			}

			for (int j = 0; j < p->macierz->mRzedow; j++) {
				PEX->A[g->elements[nrElementu]->nodesID[j] - 1][0] += p->macierz->A[j][0];
			}

		}

		Matrix* RES = new Matrix(ILOSC_WEZLOW, ILOSC_WEZLOW, 0.0);

		for (int i = 0; i < ILOSC_WEZLOW; i++)
			for (int j = 0; j < ILOSC_WEZLOW; j++) {

				RES->A[i][j] = HEX->A[i][j] + (CEX->A[i][j] / DELTA_T);
			}

		//RES->print();

		Matrix* PES = new Matrix(1, ILOSC_WEZLOW, 0.0);

		for (int i = 0; i < ILOSC_WEZLOW; i++)
			for (int j = 0; j < ILOSC_WEZLOW; j++) {
				double a = g->nodes[j]->t0;
				double b = (CEX->A[i][j] / DELTA_T);
				double c = HEX->A[i][j] + (CEX->A[i][j] / DELTA_T);
				PES->A[i][0] += g->nodes[j]->t0 * (CEX->A[i][j] / DELTA_T);
			}

		for (int i = 0; i < ILOSC_WEZLOW; i++) {
			PES->A[i][0] += PEX->A[i][0];
		}

		//PES->print();
		std::cout << "ITERACJA " << iteracja << std::endl;
		Matrix* X = gauss(RES, PES);
		X->print();

		for (int i = 0; i < ILOSC_WEZLOW; i++) {
			double a = X->A[i][0];
			g->nodes[i]->t0 = X->A[i][0];
		}

	}

	//for (int i = 0; i < ILOSC_WEZLOW; i++) {
	//	std::cout << g->nodes[i]->t0 << std::endl;
	//}


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

	return X;
}