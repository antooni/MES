#pragma once

#include "kwadratury.h"

#include <vector>

struct Element4 {
	std::vector<std::vector<double>> dKsi;
	std::vector<std::vector<double>> dEta;

	std::vector<double(*)(double)> pochodneKsi;
	std::vector<double(*)(double)> pochodneEta;

	Kwadratura* kwadratura;

	Element4();
};