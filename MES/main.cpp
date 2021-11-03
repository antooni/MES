#include <iostream>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"
#include "jakobian.h"

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
	/*Schemat s = TRZY_PUNKTOWY;
	Parametry* parametry = new Parametry(s);

	Calka* calka = new Calka(funkcja2, parametry);

	Gauss* gauss = new Gauss(calka);

	std::cout << gauss->oblicz();*/

	Element4* element = new Element4(DWU_PUNKTOWY);

	Grid* g = new Grid(0.025, 0.025, 2, 2);

	Jakobian* j = new Jakobian(0, 0, *element, *g);
	
	double r1 = j->J_inv[0] * element->dKsi[0][0] + j->J_inv[1] * element->dKsi[0][1];

	return 0;
}