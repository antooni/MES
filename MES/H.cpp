#include "H.h"
#include "jakobian.h"

constexpr auto CONDUCTIVITY = 25.0;
constexpr auto ALFA = 300.0;

H::H(int nrElementu, Element4 element, Grid grid) {
	int n_funkcji_ksztaltu = element.N_FUNKCJI;
	int m_punktow_calkowania = element.indeksyPunktow.size();

	Matrix* dNdx = new Matrix(n_funkcji_ksztaltu, m_punktow_calkowania, 0.0);
	Matrix* dNdy = new Matrix(n_funkcji_ksztaltu, m_punktow_calkowania, 0.0);

	Matrix** HPc = new Matrix * [m_punktow_calkowania]();

	for (int k = 0; k < m_punktow_calkowania; k++) {
		//tworzymy macierz HPc dla punktu calkowania
		HPc[k] = new Matrix(n_funkcji_ksztaltu, n_funkcji_ksztaltu, 0.0);
		//obliczamy jakobian dla ka¿dego punktu ca³kowania
		Jakobian* jakobian = new Jakobian(nrElementu, k, element, grid);

		//wyliczamy wartosci dNdx dNdY dla punktu calkowania (rz¹d)
		for (int j = 0; j < element.N_FUNKCJI; j++) {
			dNdx->A[k][j] = jakobian->J_inv[0] * element.dKsi[k][j] + jakobian->J_inv[1] * element.dEta[k][j];
			dNdy->A[k][j] = jakobian->J_inv[2] * element.dKsi[k][j] + jakobian->J_inv[3] * element.dEta[k][j];
		}

		//obliczamy macierz H w poszczególnych punktach ca³kowania

		for (int i = 0; i < element.N_FUNKCJI; i++) {
			for (int j = 0; j < element.N_FUNKCJI; j++) {
				double a = dNdx->A[k][i] * dNdx->A[k][j];
				double c = dNdy->A[k][i] * dNdy->A[k][j];
				//transponujemy w kolejnych komorkach
				double T_komorka = a + c;

				//pobieramy wagê
				double waga = element.punktyCalkowania[element.indeksyPunktow[k]->xIndex]->A;
				waga *= element.punktyCalkowania[element.indeksyPunktow[k]->yIndex]->A;

				//zapisujemy kolejne komorki HPc
				HPc[k]->A[i][j] = CONDUCTIVITY * (T_komorka)*jakobian->det * waga;
			}
		}
		delete jakobian;
	}

	//sumujemy macierze HPc z punktów ca³kowania
	Matrix* H = new Matrix(n_funkcji_ksztaltu, n_funkcji_ksztaltu, 0.0);

	for (int k = 0; k < m_punktow_calkowania; k++) {
		for (int i = 0; i < n_funkcji_ksztaltu; i++) {
			for (int j = 0; j < n_funkcji_ksztaltu; j++) {
				H->A[i][j] += HPc[k]->A[i][j];
			}
		}
	}
	macierz = H;

	for (int i = 0; i < m_punktow_calkowania; i++) delete HPc[i];
	delete HPc;
}

Hbc::Hbc(int nrElementu, Element4 element, Grid grid)
{
	std::vector<Walls> walls;
	//definiujemy œciany
	for (int i = 0; i < grid.elements[nrElementu]->nodesID.size(); i++) {
		int tmpNodeID = grid.elements[nrElementu]->nodesID[i] - 1;
		int tmpNodeID2;

		Node* tmpNode = grid.nodes[tmpNodeID];
		
		if (i == grid.elements[nrElementu]->nodesID.size() - 1) {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[0] - 1;
		}
		else {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[i + 1] - 1;
		}
		Node* tmpNode2 = grid.nodes[tmpNodeID2];

		if(tmpNode->isBC && tmpNode2->isBC) {
			if (i == 0) walls.push_back(BOTTOM);
			if (i == 1) walls.push_back(RIGHT);
			if (i == 2) walls.push_back(TOP);
			if (i == 3) walls.push_back(LEFT);
		}
	}

	int n = element.N.size();
	int m = element.indeksyPunktow.size();
	int p = sqrt(m);

	Matrix* HBC = new Matrix(n, n, 0.0);

	for (int nrSciany = 0; nrSciany < walls.size(); nrSciany++) {
		Matrix* nValue = new Matrix(n, p, 0.0);

		for (int i = 0; i < p; i++) {
			for (int j = 0; j < element.N_FUNKCJI; j++) {
				int wallOffset = p * walls[nrSciany];
				nValue->A[i][j] = element.N[j](element.punktyPoPowierzchni[i + wallOffset]->x, element.punktyPoPowierzchni[i + wallOffset]->y);
			}
		}
		//nValue->print();

		int tmpNodeID = grid.elements[nrElementu]->nodesID[nrSciany] - 1;
		int tmpNodeID2;

		Node* tmpNode = grid.nodes[tmpNodeID];

		if (nrSciany == 3) {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[0] - 1;
		}
		else {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[nrSciany + 1] - 1;
		}
		Node* tmpNode2 = grid.nodes[tmpNodeID2];

		int wallOffset = 2 * walls[nrSciany];
		double x1 = tmpNode->x;
		double x2 = tmpNode2->x;
		double y1 = tmpNode->y;
		double y2 = tmpNode2->y;

		double L = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
		double detJ = L / 2;

		for (int i = 0; i < element.N_FUNKCJI; i++) {
			for (int j = 0; j < element.N_FUNKCJI; j++) {
				double T_komorka = 0.0;

				for (int jj = 0; jj < p; jj++) {
					int wallOffset = p * walls[nrSciany];
					T_komorka += nValue->A[jj][i] * nValue->A[jj][j] * element.punktyPoPowierzchni[jj+wallOffset]->w;
				}

				HBC->A[i][j] += ALFA * T_komorka * detJ;
			}
		}

		delete nValue;
	}
	macierz = HBC;
}

P::P(int nrElementu, Element4 element, Grid grid)
{
	std::vector<Walls> walls;
	for (int i = 0; i < grid.elements[nrElementu]->nodesID.size(); i++) {
		int tmpNodeID = grid.elements[nrElementu]->nodesID[i] - 1;
		int tmpNodeID2;

		Node* tmpNode = grid.nodes[tmpNodeID];

		if (i == grid.elements[nrElementu]->nodesID.size() - 1) {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[0] - 1;
		}
		else {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[i + 1] - 1;
		}
		Node* tmpNode2 = grid.nodes[tmpNodeID2];

		if (tmpNode->isBC && tmpNode2->isBC) {
			if (i == 0) walls.push_back(BOTTOM);
			if (i == 1) walls.push_back(RIGHT);
			if (i == 2) walls.push_back(TOP);
			if (i == 3) walls.push_back(LEFT);
		}

	}


	int n = element.N.size();
	int m = element.indeksyPunktow.size();
	int p = sqrt(m);


	Matrix* P = new Matrix( 1, 4, 0.0);

	for (int nrSciany = 0; nrSciany < walls.size(); nrSciany++) {

		Matrix* nValue = new Matrix(n, p, 0.0);

		for (int i = 0; i < p; i++) {
			for (int j = 0; j < 4; j++) {
				int wallOffset = p * walls[nrSciany];
				nValue->A[i][j] = element.N[j](element.punktyPoPowierzchni[i + wallOffset]->x, element.punktyPoPowierzchni[i + wallOffset]->y);
			}
		}
		//nValue->print();

		int tmpNodeID = grid.elements[nrElementu]->nodesID[nrSciany] - 1;
		int tmpNodeID2;

		Node* tmpNode = grid.nodes[tmpNodeID];

		if (nrSciany == 3) {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[0] - 1;
		}
		else {
			tmpNodeID2 = grid.elements[nrElementu]->nodesID[nrSciany + 1] - 1;
		}
		Node* tmpNode2 = grid.nodes[tmpNodeID2];

		int wallOffset = 2 * walls[nrSciany];
		double x1 = tmpNode->x;
		double x2 = tmpNode2->x;
		double y1 = tmpNode->y;
		double y2 = tmpNode2->y;

		double L = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

		double detJ = L / 2;

		double* T = new double[4];

		T[BOTTOM] = 1200.0;
		T[RIGHT] = 1200.0;
		T[TOP] = 1200.0;
		T[LEFT] = 1200.0;

		for (int j = 0; j < element.N_FUNKCJI; j++) {

			double T_komorka = 0.0;

			for (int jj = 0; jj < p; jj++) {
				int wallOffset = p * walls[nrSciany];
				T_komorka += nValue->A[jj][j] * T[walls[nrSciany]] * element.punktyPoPowierzchni[jj + wallOffset]->w;
			}

			P->A[j][0] += ALFA * T_komorka * detJ ;
		}

		delete nValue;
		delete T;

	}

	macierz = P;
}
