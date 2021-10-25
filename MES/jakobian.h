#pragma once

#include "kwadratury.h"
#include "calka.h"

#include <vector>

struct Point {
	int xIndex;
	int yIndex;

	Point(int _x, int _y) {
		xIndex = _x;
		yIndex = _y;
	}
};

struct Element4 {
	std::vector<std::vector<double>> dKsi;
	std::vector<std::vector<double>> dEta;

	std::vector<double(*)(double)> pochodneKsi;
	std::vector<double(*)(double)> pochodneEta;

	std::vector<Point*> indeksyPunktow;

	Kwadratura* kwadratura;

	Element4(Schemat s);

	void wyswietlDKsi();
	void wyswietlDEta();
};