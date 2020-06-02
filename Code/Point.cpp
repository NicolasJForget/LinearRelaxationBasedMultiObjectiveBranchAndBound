#include "Point.h"

/* ==========================================================
		Constructors
 ========================================================= */

Point::Point(std::vector<double> const& z) : objVector(z), preImage(0), adjList(0), activeHyperplanes(0), discarded(false) {};

Point::Point() : objVector(0), preImage(0), adjList(0), activeHyperplanes(0), discarded(false) {};

/* ==========================================================
		Regular methods
 ========================================================= */

int Point::sizeIntersectionActiveHyperplanes(Point& u) {

	int s = 0;
	std::list<int>::iterator itThis;
	std::list<int>::iterator itu;
	std::list<int>* actHu = u.get_activeHyperplanes(); // pointer to list of active Hyperplane of u

	for (itThis = activeHyperplanes.begin(); itThis != activeHyperplanes.end(); ++itThis) {
		for (itu = actHu->begin(); itu != actHu->end(); ++itu) {
			if (*itThis == *itu) {
				++s;
			}
		}
	}

	return s;
}

bool Point::is_extreme() {

	bool extr = true;
	std::list<int>::iterator it = activeHyperplanes.begin();

	while (extr && it != activeHyperplanes.end()) {
		if (*it <= -1) { // if hyperplane is articifial
			extr = false;
		}
		++it;
	}

	return extr;
}

 /* ==========================================================
		 Getters
  ========================================================= */

double Point::get_objVector(int obj) {
	return objVector[obj];
}

bool Point::isDiscarded() {
	return discarded;
}

std::list<Point*>* Point::get_adjList() {
	return &adjList;
}

std::list<int>* Point::get_activeHyperplanes() {
	return &activeHyperplanes;
}

Point* Point::get_adress() {
	return this;
}

/* ==========================================================
		Setters
 ========================================================= */

void Point::set_objVector(int obj, double val) {
	objVector[obj] = val;
}

void Point::add_adjecentPoint(Point* adj) {
	adjList.push_back(adj);
}

void Point::add_activeHyperplane(int identifier) {
	activeHyperplanes.push_back(identifier);
}

void Point::becomesDiscarded() {
	discarded = true;
}

void Point::replaceAdjVertex(Point* oldVertex, Point* newVertex) {
	adjList.remove(oldVertex);
	adjList.push_back(newVertex);
}

void Point::updateActiveHyperplanes(Point& u, Point& v, int hyperplaneID) {

	std::list<int>* hu = u.get_activeHyperplanes();
	std::list<int>::iterator itu;
	std::list<int>* hv = v.get_activeHyperplanes();
	std::list<int>::iterator itv;

	for (itu = hu->begin(); itu != hu->end(); ++itu) {
		for (itv = hv->begin(); itv != hv->end(); ++itv) {
			if (*itu == *itv) {
				activeHyperplanes.push_back(*itu);
			}
		}
	}
	activeHyperplanes.push_back(hyperplaneID);
}