#pragma once
#include <vector>

struct PunktKwadratury {
	double x;
	double A;

	PunktKwadratury(double, double);
};

struct StopienKwadratury {
	std::vector<PunktKwadratury*> punkty;

	StopienKwadratury(std::vector<PunktKwadratury*>);
};

struct Kwadratura {
	std::vector<StopienKwadratury*> stopnie;

	Kwadratura();
};