#pragma once

#include <vector>
#include <iostream> 


struct Node
{
	double x;
	double y;
	bool isBC;

	Node(double _x, double _y, bool isBC);
};

struct Element
{
	std::vector<int> nodesID;

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

