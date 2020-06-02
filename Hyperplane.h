#pragma once
#include <vector>
#include "Point.h"
class Hyperplane
{
private:
	std::vector<double> normalVector;
	int dim; // p - 1
	double d;
	int id;

public:
	Hyperplane(std::vector<double> const& nv, std::vector<double> const& pt, int identifier);

	std::vector<Point> GenerateIntersectionTriangle(std::vector<double>& cone);
	double FindMissingCoordinate(Point& pt, int coord);
	bool Above(Point & z); // true if z is above hyperplane
	std::vector<double> EdgeIntersection(Point& u, Point& v);

	double get_d();
	int get_dim();
	double get_normalVector(int coord);
	int get_id();
};

