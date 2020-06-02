#include "LowerBoundSet.h"

/* ==========================================================
		Constructors
 ========================================================= */

LowerBoundSet::LowerBoundSet(LinearProgram linprog) : lp(linprog), antiIdealPoint(linprog.get_p()), interiorPoint(linprog.get_p()), facets(), nextId(0), extrPoints(), bensonVertexIterator(extrPoints.begin()), warmstarted(false)
{
	initialize();
}; // no warmstarting

/* ==========================================================
		Regular Methods
 ========================================================= */

int LowerBoundSet::newId() {
	++nextId;
	return nextId - 1;
}

void LowerBoundSet::initialize() {

	//antiIdealpoint init -> solve max z(x)
	for (int k = 0; k < lp.get_p(); k++) {
		antiIdealPoint[k] = 0;
	}
	
	// solve z1(x) + ... + zp(z)

	std::vector<double> normalVector(lp.get_p(), 1);
	std::vector<double> ws(lp.get_p());
	for (int i = 0; i < lp.get_p(); i++) {
		ws[i] = double(-2) - i; // to be changed, solve min z1(x) + ... + zp(z)
	}

	// add the new hyperplan to the LB set

	Hyperplane h(normalVector,ws,newId());
	facets.push_back(h);

	// update the set of points
	extrPoints.push_back(Point(antiIdealPoint));
	for (int k = 0; k < lp.get_p(); k++) {
		extrPoints.push_back(Point(antiIdealPoint));
		extrPoints.back().set_objVector(k, h.FindMissingCoordinate(extrPoints.back(), k));
		extrPoints.back().add_activeHyperplane(h.get_id());
	}

	// create adjacency lists for extreme points
	std::list<Point>::iterator it1;
	int k = 0;
	std::list<Point>::iterator it2;
	int l;
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		l = 0;
		for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
			if (k != l) {
				it1->add_adjecentPoint(it2->get_adress());
			}
			l++;
		}
		k++;
	}

	// add p artificial dominated plans (search cone)
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		for (int k = 1; k <= lp.get_p(); k++) {
			it1->add_activeHyperplane(-k);
		}
	}
}

void LowerBoundSet::print() {

	std::cout << "============= Lower Bound Set ==============\n\n";

	std::cout << "Anti-ideal point: ( ";
	for (int k = 0; k < lp.get_p() - 1; k++) {
		std::cout << antiIdealPoint[k] << " , ";
	}
	std::cout << antiIdealPoint[lp.get_p() - 1] << " )\n";

	// Hyperplane representation
	std::cout << "\nHyperplane representation:\n";
	std::list<Hyperplane>::iterator it;
	for (it = facets.begin(); it != facets.end(); ++it) {
		std::cout << "( ";
		for (int k = 0; k < lp.get_p()-1; k++) {
			std::cout << (*it).get_normalVector(k) << " , ";
		}
		std::cout << (*it).get_normalVector(lp.get_p()-1) << " ) " << " = " << (*it).get_d() << std::endl;
	}

	// Extreme points representation

	std::cout << "\n Extreme points representation:\n";
	std::list<Point>::iterator it2;
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		std::cout << " ( ";
		for (int k = 0; k < lp.get_p() - 1; k++) {
			std::cout << (*it2).get_objVector(k) << " , ";
		}
		std::cout << (*it2).get_objVector(lp.get_p()-1) << " )\n";
	}
	std::cout << "\n";
}

void LowerBoundSet::compute() {

}

void LowerBoundSet::updatePolytope(Hyperplane& H) {

	std::cout << "\nHyperplane added";
	facets.push_back(H);

	// search for infeasible points
	std::list<Point>::iterator it = extrPoints.begin();
	//std::list<Point*> discardedPts;
	while (it != extrPoints.end()) {
		if (!H.Above(*it)) { // point below H, so it is not in the polytope anymore and it should be discarded
			it->becomesDiscarded();
			std::cout << "\npoint discarded : ";
			for (int i = 0; i < lp.get_p(); i++) {
				std::cout << " " << it->get_objVector(i);
			}
			//discardedPts.push_back(it->get_adress());
		}
		else {
			std::cout << "\nnot this one...";
		}
		++it;
	}

	// search for new extreme points (comparison with edges, by adjacency lists)
	std::list<Point*>::iterator itAdj;
	std::list<Point*>* adjacentVertex;
	Point newPts;
	std::vector<Point*> allNewVertices(0);

	for (it = extrPoints.begin(); it != extrPoints.end(); ++it) { // for each extreme point
		if (it->isDiscarded()) { // if it is discarded

			std::cout << "\nnew infeasible vertex :";
			for (int i = 0; i < lp.get_p(); i++) {
				std::cout << " " << it->get_objVector(i);
			}

			adjacentVertex = it->get_adjList(); // we look at its adjacency list
			for (itAdj = adjacentVertex->begin(); itAdj != adjacentVertex->end(); ++itAdj) {
				
				std::cout << "\ntest with :";
				for (int i = 0; i < lp.get_p(); i++) {
					std::cout << " " << (*itAdj)->get_objVector(i);
				}

				if (!(*itAdj)->isDiscarded()) { // for each non-discarded adjacent vertex, we compute the new point

					std::cout << "\n   -> Non-discarded adjecent vertex :";
					for (int i = 0; i < lp.get_p(); i++) {
						std::cout << " " << (*itAdj)->get_objVector(i);
					}

					newPts = Point( H.EdgeIntersection(*it,**itAdj) ); // new pts as intersection of H and edge
					(*itAdj)->replaceAdjVertex(it->get_adress(), &newPts); // update adj vertex for feasible one
					newPts.add_adjecentPoint((*itAdj)->get_adress()); // init adj vertex for new one
					newPts.updateActiveHyperplanes(*it,**itAdj,H.get_id()); // compute active hyperplanes
					extrPoints.push_back(newPts); //new pts added to list of extreme points
					allNewVertices.push_back(&newPts); // remember the new vertex
				}
			}
		}
	}

	// Update the adjacency lists for the new vertices
	int s = allNewVertices.size();
	for (int i = 0; i < s - 1; i++) {
		for (int j = i + 1; j < s; j++) {
			if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) == H.get_dim()) {
				allNewVertices[i]->add_adjecentPoint(allNewVertices[j]);
				allNewVertices[j]->add_adjecentPoint(allNewVertices[i]);
			}
		}
	}

	// delete non-feasible vertices
	std::list<Point>::iterator itdelete = extrPoints.begin();
	std::list<Point>::iterator current;
	while (itdelete != extrPoints.end()) {
		current = itdelete;
		++itdelete;
		if (current->isDiscarded()) {
			extrPoints.erase(current);
		}
	}

	std::cout << "\n\n";
}

/* ==========================================================
		Getters
 ========================================================= */

 /* ==========================================================
		 Setters
  ========================================================= */