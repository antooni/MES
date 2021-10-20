#include "jakobian.h"

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


Element4::Element4()
{
	dKsi = {
		{},{},{},{}
	};

	dEta = {
		{},{},{},{}
	};

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

	std::vector<PunktKwadratury*> punkty = kwadratura->stopnie[0]->punkty;


	for (int i = 0; i < pochodneKsi.size(); i++) {
		dKsi[0].push_back(pochodneKsi[i](punkty[0]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dKsi[1].push_back(pochodneKsi[i](punkty[0]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dKsi[2].push_back(pochodneKsi[i](punkty[1]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dKsi[3].push_back(pochodneKsi[i](punkty[1]->x));
	}






	for (int i = 0; i < pochodneKsi.size(); i++) {
		dEta[0].push_back(pochodneEta[i](punkty[0]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dEta[1].push_back(pochodneEta[i](punkty[1]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dEta[2].push_back(pochodneEta[i](punkty[1]->x));
	}

	for (int i = 0; i < pochodneKsi.size(); i++) {
		dEta[3].push_back(pochodneEta[i](punkty[0]->x));
	}
	
}
