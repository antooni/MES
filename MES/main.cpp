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

double fun(double x) {return 5 * x * x + 3 * x + 6;}

double fun2(double x, double y) {return 5 * x * x * y * y + 3 * x * y + 6;}

Funkcja* funkcja = new Funkcja(fun, "f(x) = 5x^2 + 3x + 6");
Funkcja* funkcja2 = new Funkcja(fun2, "f(x) = 5x^2y^2 + 3xy + 6");

bool gauss(int n, double** AB, double* X);


int main()
{
	Element4* element = new Element4(DWU_PUNKTOWY);

	Grid* g = new Grid(0.1, 0.1, 4, 4);

	g->homework();

	const int ILOSC_WEZLOW = g->nodes.size();
	
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

	Matrix* AB = new Matrix(ILOSC_WEZLOW + 1, ILOSC_WEZLOW, 0.0);

	for (int i = 0; i < AB->nKolumn; i++) {
		for (int j = 0; j < AB->mRzedow; j++) {
			if (i == 16) AB->A[j][i] = PEX->A[j][0];
			else AB->A[i][j] = HEX->A[i][j];
		}
	}
	AB->print();

	double* X = new double[ILOSC_WEZLOW];
	if (gauss(ILOSC_WEZLOW, AB->A, X))
	{
		for (int i = 0; i < ILOSC_WEZLOW; i++)
			std::cout << "x" << i + 1 << " = "  << X[i]
			<< std::endl;
	}
	else
		std::cout << "DZIELNIK ZERO\n";

	CEX->print();

	return 0;
}

bool gauss(int n, double** AB, double* X)
{
	int i, j, k;
	double m, s;

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}