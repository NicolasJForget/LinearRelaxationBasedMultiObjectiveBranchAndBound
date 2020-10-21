#include <stdlib.h>
#include "LinearProgram.h"
#include "Hyperplane.h"
#include <list>
#include <vector>
#include "Point.h"

#pragma once
class LowerBoundSet
{
private:
	LinearProgram lp;
	CplexModel* weightedSum;
	CplexModel* feasibilityModel;
	CplexModel* dualBenson;
	CplexModel* bestValidPoint;

	std::vector <double> antiIdealPoint;
	std::vector <double> interiorPoint;

	std::list <Hyperplane*> facets;
	std::vector <Hyperplane*> boundingBox;
	int nextHyperplaneId;
	int nextPointId;
	std::list <Point*> extrPoints;
	std::list<Point*>::iterator bensonVertexIterator;

	bool warmstarted;
	int nbIterations;

public:

	LowerBoundSet(MathematicalFormulation linprog, CplexModel* ws, CplexModel* feas, CplexModel* db, CplexModel* bvp); // Constructor, no warmstarting

	int newHyperplaneId(); // Give a new id to an hyperplane
	int newPointId(); // Give a new id to a point
	void initialize(MathematicalFormulation& linprog); // Simplex initialization for no warmstarting
	void print(); // print the lower bound set
	void compute(MathematicalFormulation& LP); // compute the lower bound set corresponding to lp
	void updatePolytope(Hyperplane* H); // update the polytope of the LB set with the new hyperplane
	void classifyExtremePoints();
};

