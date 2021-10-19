#pragma once

#include "kwadratury.h"
#include <vector>
#include <cmath>

Kwadratura::Kwadratura()
{
	std::vector<PunktKwadratury*> punkty1 = { 
		new PunktKwadratury(-1 / sqrt(3.0),1), 
		new PunktKwadratury(1 / sqrt(3.0),1), 
	};
	std::vector<PunktKwadratury*> punkty2 = {
		new PunktKwadratury(-sqrt((3.0/5.0)),(5.0/9.0)),
		new PunktKwadratury(0,(8.0/9.0)),
		new PunktKwadratury(sqrt(3.0/5.0),(5.0/9.0)),
	};
	std::vector<PunktKwadratury*> punkty3 = {
		new PunktKwadratury(-0.861136,0.347855),
		new PunktKwadratury(-0.339981,0.652145),
		new PunktKwadratury(0.339981,0.652145),
		new PunktKwadratury(0.861136,0.347855),

	};
	stopnie = {
		new StopienKwadratury(punkty1),
		new StopienKwadratury(punkty2),
		new StopienKwadratury(punkty3)
	};
}

StopienKwadratury::StopienKwadratury(std::vector<PunktKwadratury*> _punkty)
{
	punkty = _punkty;
}

PunktKwadratury::PunktKwadratury(double _x, double _A)
{
	x = _x;
	A = _A;
}