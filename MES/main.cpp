#include <iostream>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"
#include "jakobian.h"
#include "matrix.h"
#include "H.h"

const double eps = 1e-12;

double fun(double x) {
	return 5 * x * x + 3 * x + 6;
}

double fun2(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

Funkcja* funkcja = new Funkcja(fun, "f(x) = 5x^2 + 3x + 6");
Funkcja* funkcja2 = new Funkcja(fun2, "f(x) = 5x^2y^2 + 3xy + 6");

bool gauss(int n, double** AB, double* X)
{

	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

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

	// wyliczanie niewiadomych

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


int main()
{
	//Schemat s = TRZY_PUNKTOWY;
	//Parametry* parametry = new Parametry(s);

	//Calka* calka = new Calka(funkcja2, parametry);

	//Gauss* gauss = new Gauss(calka);

	//std::cout << gauss->oblicz();

	Element4* element = new Element4(DWU_PUNKTOWY);

	//Grid* g = new Grid(0.2, 0.1, 5, 4);

	//Grid* g = new Grid(0.025, 0.025, 2, 2);
	Grid* g = new Grid(0.1, 0.1, 4, 4);
	
	
	//Grid* g = new Grid(0.025, 0.025, 6, 2);



	g->homework();

	Hbc* hbc = new Hbc(0, *element, *g);


	Matrix* HEX = new Matrix(16, 16, 0.0);
	Matrix* PEX = new Matrix(1, 16, 0.0);

	for (int nrElementu = 0; nrElementu < g->elements.size(); nrElementu++) {
		H* h = new H(nrElementu, *element, *g);
		Hbc* hbc = new Hbc(nrElementu, *element, *g);
		P* p = new P(nrElementu, *element, *g);

		Matrix res = h->macierz->add(hbc->macierz);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				HEX->A[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += res.A[i][j];
			}
		}

		for (int j = 0; j < 4; j++) {
			// TODO  agregacja P
			PEX->A[g->elements[nrElementu]->nodesID[j] - 1][0] += p->macierz->A[j][0];
		}

	}

	//HEX->print();
	//PEX->print();

	int n = 16;

	Matrix* AB = new Matrix(17, 16, 0.0);


	for (int i = 0; i < 17; i++) {
		for (int j = 0; j < 16; j++) {
			if (i == 16) AB->A[j][i] = PEX->A[j][0];
			else AB->A[i][j] = HEX->A[i][j];
		}
	}
	AB->print();

	double* X = new double[n];

	if (gauss(n, AB->A, X))
	{
		for (int i = 0; i < n; i++)
			std::cout << "x" << i + 1 << " = "  << X[i]
			<< std::endl;
	}
	else
		std::cout << "DZIELNIK ZERO\n";

	//HEX->print();

	//AB->print();






	return 0;
}

