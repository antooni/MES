#pragma once

#include <vector>
#include <iostream> 


struct Node
{
	double x;
	double y;

	Node(double _x, double _y);
};

struct Element
{
	std::vector<int> id;

	Element(std::vector<int> _id);
};

struct Grid
{
	std::vector<Node*> nodes;
	std::vector<Element*> elements;

	double H;
	double B;
	int nH;
	int nB;
	int nN;
	int nE;

	Grid(double _H, double _B, int _nH, int _nN);

	void homework();


};

