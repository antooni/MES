#pragma once

#include "kwadratury.h"
#include "calka.h"
#include "grid.h"

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

struct Element4 {
	std::vector<std::vector<double>> dKsi;
	std::vector<std::vector<double>> dEta;

	std::vector<double(*)(double)> pochodneKsi;
	std::vector<double(*)(double)> pochodneEta;

	std::vector<IndexPoint*> indeksyPunktow;

	Kwadratura* kwadratura;

	Element4(Schemat s);

	void wyswietlDKsi();
	void wyswietlDEta();
};



struct Jakobian {
	std::vector<double> J;
	std::vector<double> J_inv;

	Jakobian(int i, int j, Element4 element, Grid grid){
		double tmp = 0;

		for (int k = 0; k < grid.elements[i]->nodesID.size(); k++) {
			//double a = element.dKsi[j][k];
			//double b = (grid.nodes[grid.elements[i]->nodesID[k] - 1]->x);

			tmp += element.dKsi[j][k] * (grid.nodes[grid.elements[i]->nodesID[k]-1.0]->x);
		}
		J.push_back(tmp);

		tmp = 0;
		for (int k = 0; k < grid.elements[i]->nodesID.size(); k++) {
			tmp += element.dKsi[j][k] * (grid.nodes[grid.elements[i]->nodesID[k] - 1.0]->y);
		}
		J.push_back(tmp);

		tmp = 0;
		for (int k = 0; k < grid.elements[i]->nodesID.size(); k++) {
			tmp += element.dEta[j][k] * (grid.nodes[grid.elements[i]->nodesID[k] - 1.0]->x);
		}
		J.push_back(tmp);

		tmp = 0;
		for (int k = 0; k < grid.elements[i]->nodesID.size(); k++) {
			tmp += element.dEta[j][k] * (grid.nodes[grid.elements[i]->nodesID[k] - 1.0]->y);
		}
		J.push_back(tmp);

		double detJ = (J[0] * J[3]) - (J[2] * J[3]);
		detJ = 1.0 / detJ;

		J_inv.push_back(detJ * J[3]);
		J_inv.push_back(-detJ * (J[1]));
		J_inv.push_back(-detJ * (J[2]));
		J_inv.push_back(detJ * J[0]);
	}

};