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

struct WPoint {
	double x;
	double y;
	double w;

	WPoint(double _x, double _y, double _w) {
		x = _x;
		y = _y;
		w = _w;
	}
};

struct Element4 {
	std::vector<std::vector<double>> dKsi;
	std::vector<std::vector<double>> dEta;

	std::vector<double(*)(double)> pochodneKsi;
	std::vector<double(*)(double)> pochodneEta;
	std::vector<double(*)(double, double)> N;

	std::vector<PunktKwadratury*> punktyCalkowania;
	std::vector<IndexPoint*> indeksyPunktow;
	std::vector<WPoint*> punktyPoPowierzchni;

	const int WYMIAR = 2;
	const int N_FUNKCJI = 4;

	Kwadratura* kwadratura;

	Schemat s;

	Element4(Schemat s);

	void wyswietlDKsi();
	void wyswietlDEta();
};



struct Jakobian {
	double* J;
	double* J_inv;

	double det;

	Jakobian(int i, int j, Element4 element, Grid grid);
	~Jakobian() {
		delete J;
		delete J_inv;
	}
};

