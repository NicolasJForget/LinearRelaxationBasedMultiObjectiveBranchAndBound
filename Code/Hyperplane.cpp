#include "Hyperplane.h"

/* ==========================================================
		Constructors
 ========================================================= */

Hyperplane::Hyperplane(std::vector<double> const& nv, std::vector<double> const& pt, int identifier) : normalVector(nv), dim(nv.size() - 1), id(identifier) {

	d = 0;
	for (int i = 0; i < dim + 1; i++) {
		d += nv[i] * pt[i];
	}
}

/* ==========================================================
		Regular Methods
 ========================================================= */

double Hyperplane::FindMissingCoordinate(Point& pt, int coord) {

	//objVector[coord] = h.get_d();
	double val = d;
	for (int k = 0; k < dim + 1; k++) {
		if (k != coord) {
			val -= normalVector[k] * pt.get_objVector(k);
		}
	}
	val /= normalVector[coord];

	return val;
}

std::vector<Point> Hyperplane::GenerateIntersectionTriangle(std::vector<double>& cone) {

	std::vector<Point> triangle(dim + 1);
	for (int k = 0; k < dim + 1; k++) {
		triangle[k] = Point(cone);
		triangle[k].set_objVector(k,FindMissingCoordinate(triangle[k], k));
	}

	return triangle;
} // useless ?

bool Hyperplane::Above(Point & z) {

	double val = 0;
	for (int k = 0; k < dim + 1; k++) {
		val += normalVector[k] * z.get_objVector(k);
	}

	return val >= d;
}

std::vector<double> Hyperplane::EdgeIntersection(Point& u, Point& v) {

	// compute the value of lambda, to find the point on the edge
	double denominator = 0;
	double lambda = d;
	for (int l = 0; l < dim + 1; l++) {
		lambda -= normalVector[l] * v.get_objVector(l);
		denominator += normalVector[l] * (u.get_objVector(l) - v.get_objVector(l));
	}
	lambda /= denominator;

	// compute the actual point
	std::vector<double> intersection(dim + 1);
	for (int k = 0; k < dim + 1; k++) {
		intersection[k] = lambda * u.get_objVector(k) + (1 - lambda) * v.get_objVector(k);
	}

	return intersection;
}

/* ==========================================================
		Getters
 ========================================================= */

double Hyperplane::get_d() {
	return d;
}

int Hyperplane::get_dim() {
	return dim;
}

double Hyperplane::get_normalVector(int coord) {
	return normalVector[coord];
}

int Hyperplane::get_id() {
	return id;
}

/* ==========================================================
		Setters
 ========================================================= */