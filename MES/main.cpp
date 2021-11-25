#include <iostream>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"
#include "jakobian.h"
#include "matrix.h"
#include "H.h"


double fun(double x) {
	return 5 * x * x + 3 * x + 6;
}

double fun2(double x, double y) {
	return 5 * x * x * y * y + 3 * x * y + 6;
}

Funkcja* funkcja = new Funkcja(fun, "f(x) = 5x^2 + 3x + 6");
Funkcja* funkcja2 = new Funkcja(fun2, "f(x) = 5x^2y^2 + 3xy + 6");

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

	for (int nrElementu = 0; nrElementu < g->elements.size(); nrElementu++) {
		H* h = new H(nrElementu, *element, *g);
		Hbc* hbc = new Hbc(nrElementu, *element, *g);

		Matrix res = h->macierz->add(hbc->macierz);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				HEX->a[g->elements[nrElementu]->nodesID[i] - 1][g->elements[nrElementu]->nodesID[j] - 1] += h->macierz->a[i][j];
			}
		}

	}

	HEX->print();





	return 0;
}