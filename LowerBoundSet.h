#include "Hyperplane.h"
#include <list>
#include <vector>
#include "Point.h"
#include "LinearProgram.h"

#pragma once
class LowerBoundSet
{
private:
	LinearProgram lp;
	std::vector <double> antiIdealPoint;
	std::vector <double> interiorPoint;
	std::list <Hyperplane> facets;
	int nextId;
	std::list <Point> extrPoints;
	std::list<Point>::iterator bensonVertexIterator;
	bool warmstarted;

public:

	LowerBoundSet(LinearProgram linprog); // Constructor, no warmstarting

	int newId(); // Give a new id to an hyperplane
	void initialize(); // Simplex initialization for no warmstarting
	void print(); // print the lower bound set
	void compute(); // compute the lower bound set corresponding to lp
	void updatePolytope(Hyperplane& H); // update the polytope of the LB set with the new hyperplane
};

