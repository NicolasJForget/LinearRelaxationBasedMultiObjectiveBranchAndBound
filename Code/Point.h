#pragma once
#include <vector>
#include <list>

class Point
{
private:
	std::vector<double> objVector;
	std::vector<double> preImage;
	std::list<Point*> adjList;
	std::list<int> activeHyperplanes;
	bool discarded; //true if the point becomes discarded by a new hyperplane added to LB set. All should be false after computation of LB set is finished.

public:
	Point(std::vector<double> const& z);
	Point();

	int sizeIntersectionActiveHyperplanes(Point& u);
	bool is_extreme();

	double get_objVector(int obj);
	bool isDiscarded();
	std::list<Point*>* get_adjList();
	std::list<int>* get_activeHyperplanes();
	Point* get_adress();

	void set_objVector(int obj, double val);
	void add_adjecentPoint(Point* adj);
	void add_activeHyperplane(int identifier);
	void becomesDiscarded();
	void replaceAdjVertex(Point* oldVertex, Point* newVertex); // can be opti in cpp if we can change value of an element in a list
	void updateActiveHyperplanes(Point& u, Point& v, int hyperplaneID);
};
