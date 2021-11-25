#include "H.h"
#include "jakobian.h"

constexpr auto CONDUCTIVITY = 25;


H::H(int nrElementu, Element4 element, Grid grid) {
	int n = element.dKsi.size();
	int m = element.indeksyPunktow.size();

	Matrix* dNdx = new Matrix(n, m, 0.0);
	Matrix* dNdy = new Matrix(n, m, 0.0);


	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			Jakobian* jakobian = new Jakobian(nrElementu, j, element, grid);

			dNdx->a[i][j] = jakobian->J_inv[0] * element.dKsi[i][j] + jakobian->J_inv[1] * element.dEta[i][j];
			dNdy->a[i][j] = jakobian->J_inv[2] * element.dKsi[i][j] + jakobian->J_inv[3] * element.dEta[i][j];
		}
	}

	Matrix** HH = new Matrix* [m]();

	for (int k = 0; k < m; k++) {
		HH[k] = new Matrix(n, m, 0.0);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				Jakobian* jakobian = new Jakobian(nrElementu, j, element, grid);

				double a = dNdx->a[k][i] * dNdx->a[k][j];

				double c = dNdy->a[k][i] * dNdy->a[k][j];

				HH[k]->a[i][j] = CONDUCTIVITY * (a + c) * jakobian->det;
			}
		}
	}

	Matrix* H = new Matrix(n, m, 0.0);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < m; k++) H->a[i][j] += HH[k]->a[i][j];
		}
	}


	std::cout << "H : " << std::endl;
	H->print();
	macierz = H;

}

Hbc::Hbc(int nrElementu, Element4 element, Grid grid)
{

	std::vector<Point*> punktyPoPowierzchni = {
		new Point(1 / sqrt(3.0), -1),
		new Point(-1 / sqrt(3.0), -1),

		new Point(1, -1 / sqrt(3.0)),
		new Point(1, 1 / sqrt(3.0)),

		new Point(1 / sqrt(3.0),1),
		new Point(-1 / sqrt(3.0),1),

		new Point(-1,1 / sqrt(3.0)),
		new Point(-1,-1 / sqrt(3.0)),
	};

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

		//works only for 2d
		if(tmpNode->isBC && tmpNode2->isBC) {
			if (i == 0) walls.push_back(BOTTOM);
			if (i == 1) walls.push_back(RIGHT);
			if (i == 2) walls.push_back(TOP);
			if (i == 3) walls.push_back(LEFT);
		}
	}


	int n = element.N.size();
	int m = element.indeksyPunktow.size();

	Matrix* HBC = new Matrix(n, m, 0.0);

	for (int nrSciany = 0; nrSciany < walls.size(); nrSciany++) {

		Matrix* nValue = new Matrix(n, 2, 0.0);

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 4; j++) {
				int wallOffset = 2 * walls[nrSciany];
				nValue->a[i][j] = element.N[j](punktyPoPowierzchni[i + wallOffset]->x, punktyPoPowierzchni[i + wallOffset]->y);
			}
		}
		nValue->print();

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

		//Matrix* czerwony = new Matrix(4, 4, 0.0);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				double a = nValue->a[0][i] * nValue->a[0][j];

				double c = nValue->a[1][i] * nValue->a[1][j];

				//czerwony->a[i][j] = (25 * a + 25 * c) * detJ;
				HBC->a[i][j] += (25 * a + 25 * c) * detJ;
			}
		}
		//czerwony->print();
	}
	macierz = HBC;
	macierz->print();
}
