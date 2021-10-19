#pragma once
#include <string>

enum class Przestrzen {
	JEDEN_D,
	DWA_D
};

enum Schemat {
	DWU_PUNKTOWY,
	TRZY_PUNKTOWY,
	CZTERO_PUNKTOWY,
};

struct Funkcja {
	std::string wzor;
	Przestrzen przestrzen;

	double (*fun_x)(double);
	double (*fun_xy)(double, double);

	double oblicz(double x);
	double oblicz(double x, double y);

	Funkcja(double (*_fun_x)(double), std::string _wzor) {
		fun_x = _fun_x;
		wzor = _wzor;
		przestrzen = Przestrzen::JEDEN_D;
	}

	Funkcja(double (*_fun_xy)(double, double), std::string _wzor) {
		fun_xy = _fun_xy;
		wzor = _wzor;
		przestrzen = Przestrzen::DWA_D;
	}

};

struct Parametry {
	Schemat schemat;

	Parametry(Schemat _s) {
		schemat = _s;
	}
};

struct Calka
{
	Funkcja* funkcja;
	Parametry* parametry;

	Calka(Funkcja* f, Parametry* p) {
		funkcja = f;
		parametry = p;
	}

};


