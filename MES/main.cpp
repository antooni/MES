#include <iostream>

#include "grid.h"
#include "kwadratury.h"
#include "calka.h"
#include "gauss.h"

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
	Schemat s = TRZY_PUNKTOWY;
	Parametry* parametry = new Parametry(s);

	Calka* calka = new Calka(funkcja2, parametry);

	Gauss* gauss = new Gauss(calka);

	std::cout << gauss->oblicz();
	
	return 0;
}