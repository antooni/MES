#include <iostream>

#include "gauss.h"
#include "kwadratury.h"
#include "calka.h"

double Gauss::oblicz()
{
	std::vector <PunktKwadratury*> punkty = kwadratura->stopnie[calka->parametry->schemat]->punkty;
	if (punkty.size() == 0) {
		std::cout << "Problem ze scheamtem calkowania, kwadratura nie istnieje" << std::endl;
		return -1.0;
	}

	double wynik = 0.0;

	switch (calka->funkcja->przestrzen)
	{
	case Przestrzen::JEDEN_D:
		for (int i = 0; i < punkty.size(); i++) {
			wynik += calka->funkcja->oblicz(punkty[i]->x) * (punkty[i]->A);
		}
		break;

	case Przestrzen::DWA_D:
		for (int i = 0; i < punkty.size(); i++) {
			for (int j = 0; j < punkty.size(); j++) {
				wynik += calka->funkcja->oblicz(punkty[i]->x, punkty[j]->x) * (punkty[i]->A) * (punkty[j]->A);
			}
		}
		break;

	default:
		std::cout << "ERR WITH PARAMS - dimensions" << std::endl;
		break;
	}

	return wynik;
}
