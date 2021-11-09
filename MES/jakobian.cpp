#include "jakobian.h"
#include "calka.h"

#include <iostream>

double pKsi_1(double eta) {
	return -0.25 * (1.0 - eta);
}

double pKsi_2(double eta) {
	return 0.25 * (1.0 - eta);
}

double pKsi_3(double eta) {
	return 0.25 * (1.0 + eta);
}

double pKsi_4(double eta) {
	return -0.25 * (1.0 + eta);
}

double pEta_1(double ksi) {
	return -0.25 * (1.0 - ksi);
}

double pEta_2(double ksi) {
	return -0.25 * (1.0 + ksi);
}

double pEta_3(double ksi) {
	return 0.25 * (1.0 + ksi);
}

double pEta_4(double ksi) {
	return 0.25 * (1.0 - ksi);
}


void Pochodne::wyswietlDKsi()
{
	std::cout << "dKsi" << std::endl;
	for (int i = 0; i < dKsi.size(); i++) {
		for (int j = 0; j < dKsi[i].size(); j++) {
			std::cout << dKsi[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Pochodne::wyswietlDEta()
{
	std::cout << "dEta" << std::endl;
	for (int i = 0; i < dEta.size(); i++) {
		for (int j = 0; j < dEta[i].size(); j++) {
			std::cout << dEta[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

Pochodne::Pochodne(Schemat s)
{
	pochodneKsi = {
		pKsi_1,
		pKsi_2,
		pKsi_3,
		pKsi_4,
	};

	pochodneEta = {
		pEta_1,
		pEta_2,
		pEta_3,
		pEta_4,
	};

	kwadratura = new Kwadratura();

	std::vector<PunktKwadratury*> punktyCalkowania;

	if (s == DWU_PUNKTOWY) {
		indeksyPunktow = {
			new IndexPoint(0,0),
			new IndexPoint(1,0),
			new IndexPoint(1,1),
			new IndexPoint(0,1),
		};

		punktyCalkowania = kwadratura->stopnie[0]->punkty;
	}
	else if (s == TRZY_PUNKTOWY) {
		indeksyPunktow = {
			new IndexPoint(0,0),
			new IndexPoint(0,1),
			new IndexPoint(0,2),

			new IndexPoint(1,0),
			new IndexPoint(1,1),
			new IndexPoint(1,2),

			new IndexPoint(2,0),
			new IndexPoint(2,1),
			new IndexPoint(2,2),
		};

		punktyCalkowania = kwadratura->stopnie[1]->punkty;
	}


	for (int i = 0; i < indeksyPunktow.size(); i++) {
		std::vector<double> tempDKsi;
		std::vector<double> tempDEta;

		for (int j = 0; j < pochodneKsi.size(); j++) {
			tempDKsi.push_back(pochodneKsi[j](punktyCalkowania[indeksyPunktow[i]->yIndex]->x));
			tempDEta.push_back(pochodneEta[j](punktyCalkowania[indeksyPunktow[i]->xIndex]->x));
		}

		dKsi.push_back(tempDKsi);
		dEta.push_back(tempDEta);
	}
}

Jakobian::Jakobian(int i, int j, Pochodne element, Grid grid) {
	double tmp = 0;

	for (int k = 0; k < grid.elements[i]->nodesID.size(); k++) {
		//double a = element.dKsi[j][k];
		//double b = (grid.nodes[grid.elements[i]->nodesID[k] - 1]->x);

		tmp += element.dKsi[j][k] * (grid.nodes[grid.elements[i]->nodesID[k] - 1.0]->x);
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
	det = detJ;
	detJ = 1.0 / detJ;

	J_inv.push_back(detJ * J[3]);
	J_inv.push_back(-detJ * (J[1]));
	J_inv.push_back(-detJ * (J[2]));
	J_inv.push_back(detJ * J[0]);
}

