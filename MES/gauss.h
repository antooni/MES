#pragma once

#include "calka.h"
#include "kwadratury.h"


struct Gauss {
	Calka* calka;
	Kwadratura* kwadratura;

	Gauss(Calka* c) {
		calka = c;
		kwadratura = new Kwadratura();
	}

	double oblicz();



};
