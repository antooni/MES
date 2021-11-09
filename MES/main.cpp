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

	Pochodne* element = new Pochodne(DWU_PUNKTOWY);

	Grid* g = new Grid(0.2, 0.1, 5, 4);
	//Grid* g = new Grid(0.025, 0.025, 2, 2);

	g->homework();

	for (int i = 0; i < g->elements.size(); i++) {
		if (i == 8) {
			std::cout << "";
		}
		Jakobian* j = new Jakobian(0, 0, *element, *g);
		H* h = new H(*j, *element, *g);
	}





	return 0;
}