#include "grid.h"

Grid::Grid(double _H, double _B, int _nH, int _nB)
{
	H = _H;
	B = _B;
	nH = _nH;
	nB = _nB;
	nN = nH * nB;
	nE = (nB - 1) * (nH - 1);

	double deltaY = H / ((double)nH - 1.0);
	double deltaX = B / ((double)nB - 1.0);

	for (double i = 0.0; i < nB; i++) {
		for (double ii = 0.0; ii < nH; ii++) {
			nodes.push_back(new Node(i * deltaX, ii * deltaY));
		}
	}

	for (int i = 1; i <= nE; i++) {
		int offset = i / nH;

		std::vector<int> ids = { i + offset,i + nH + offset,i + nH + 1 + offset, i + 1 + offset };

		elements.push_back(new Element(ids));
	}
}

void Grid::homework()
{
	std::cout << "Liczba wezlow nN: " << nN << std::endl;
	std::cout << "Liczba elementow nE: " << nE << std::endl << std::endl;

	std::cout << "=== WEZLY === (x,y)" << std::endl;
	for (int i = 0; i < nN; i++) {
		std::cout << i+1 << ": (" << nodes[i]->x << " , " << nodes[i]->y <<")"<<std::endl;
	}
	std::cout << std::endl;

	std::cout << "=== ELEMENTY === [ID1,ID2,ID3,ID4] " << std::endl;
	for (int i = 0; i  < nE; i++) {
		std::cout << i+1 << ": [" << elements[i]->id[0] << " , " << elements[i]->id[1] << " , "  << elements[i]->id[2] << " , " << elements[i]->id[3] << "]" << std::endl;
	}
	std::cout << std::endl;

}

Node::Node(double _x, double _y) {
	x = _x;
	y = _y;
}

Element::Element(std::vector<int> _id) {
	id = _id;
}



//Grid::Grid() {
//
//
//	double deltaY = H / ((double)nH - 1.0);
//	double deltaX = B / ((double)nB - 1.0);
//
//	for (double i = 0.0; i < nB; i++) {
//		for (double ii = 0.0; ii < nH; ii++) {
//			nodes.push_back(new Node(i * deltaX, ii * deltaY));
//		}
//	}
//
//	for (int i = 1; i <= nE; i++) {
//		int offset = i / nH;
//
//		std::vector<int> ids = { i + offset,i + nH + offset,i + nH + 1 + offset, i + 1 + offset };
//
//		elements.push_back(new Element(ids));
//	}
//}