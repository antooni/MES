#pragma once

#include "kwadratury.h"
#include "calka.h"
#include "grid.h"
#include "matrix.h"

#include <vector>

struct IndexPoint {
	int xIndex;
	int yIndex;

	IndexPoint(int _x, int _y) {
		xIndex = _x;
		yIndex = _y;
	}
};

struct Point {
	double x;
	double y;

	Point(double _x, double _y) {
		x = _x;
		y = _y;
	}
};

struct Pochodne {
	std::vector<std::vector<double>> dKsi;
	std::vector<std::vector<double>> dEta;

	std::vector<double(*)(double)> pochodneKsi;
	std::vector<double(*)(double)> pochodneEta;

	std::vector<IndexPoint*> indeksyPunktow;

	Kwadratura* kwadratura;

	Pochodne(Schemat s);

	void wyswietlDKsi();
	void wyswietlDEta();
};



struct Jakobian {
	std::vector<double> J;
	std::vector<double> J_inv;

	double det;

	Jakobian(int i, int j, Pochodne element, Grid grid);
};

