#include "H.h"
#include "jakobian.h"

H::H(Jakobian jakobian, Pochodne element, Grid grid) {

	int n = element.dKsi.size();
	int m = element.indeksyPunktow.size();

	Matrix* dNdx = new Matrix(n, m, 0.0);
	Matrix* dNdy = new Matrix(n, m, 0.0);


	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			dNdx->a[i][j] = jakobian.J_inv[0] * element.dKsi[i][j] + jakobian.J_inv[1] * element.dEta[i][j];
			dNdy->a[i][j] = jakobian.J_inv[2] * element.dKsi[i][j] + jakobian.J_inv[3] * element.dEta[i][j];
		}
	}

	Matrix** HH = new Matrix* [m]();

	for (int k = 0; k < m; k++) {
		HH[k] = new Matrix(n, m, 0.0);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				double a = dNdx->a[k][i] * dNdx->a[k][j];

				double c = dNdy->a[k][i] * dNdy->a[k][j];

				HH[k]->a[i][j] = 30 * (a + c) * jakobian.det;
			}
		}
	}

	Matrix* H = new Matrix(n, m, 0.0);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < m; k++) H->a[i][j] += HH[k]->a[i][j];
		}
	}

	H->print();
	macierz = H;

}