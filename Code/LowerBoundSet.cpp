#include "LowerBoundSet.h"

/* ==========================================================
		Constructors
 ========================================================= */

LowerBoundSet::LowerBoundSet(MathematicalFormulation linprog, CplexModel* ws, CplexModel* feas, CplexModel* db, CplexModel* bvp) : lp(linprog), weightedSum(ws), feasibilityModel(feas), dualBenson(db), bestValidPoint(bvp), facets(), boundingBox(linprog.get_p()), nextHyperplaneId(0), nextPointId(0), extrPoints(), bensonVertexIterator(extrPoints.begin()), warmstarted(false), nbIterations(0)
{
	// init anti-ideal & interior point
	antiIdealPoint = std::vector<double>(linprog.get_p());
	interiorPoint = std::vector<double>(linprog.get_p());
	for (int k = 0; k < linprog.get_p(); k++) {
		for (int i = 0; i < linprog.get_n(); i++) {
			if (linprog.get_objective(k, i) >= 0) {
				antiIdealPoint[k] += linprog.get_objective(k, i);
			}
		}
		antiIdealPoint[k] += 2;
		interiorPoint[k] = antiIdealPoint[k] - 1; // arbitrary, to be modified ? + do we actually need interior pts ?
	}
	initialize(linprog);

}; // no warmstarting

/* ==========================================================
		Regular Methods
 ========================================================= */

int LowerBoundSet::newHyperplaneId() {
	++nextHyperplaneId;
	return nextHyperplaneId - 1;
}

int LowerBoundSet::newPointId() {
	++nextPointId;
	return nextPointId - 1;
}

void LowerBoundSet::initialize(MathematicalFormulation& linprog) {
	
	// solve z1(x) + ... + zp(z)

	std::vector<double> normalVector(lp.get_p(), 1);
	double ws = weightedSum->getWeightedSumValue(linprog,normalVector); // ws value -> does it defines rhs of this hyperplane ?
	

	// add the new hyperplan to the LB set

	Hyperplane* h = new Hyperplane(normalVector,ws);
	facets.push_back(h);

	// update the set of points
	Point* newPts = new Point(antiIdealPoint);
	extrPoints.push_back(newPts);
	//h->addVertex();
	for (int k = 0; k < lp.get_p(); k++) {
		extrPoints.push_back(new Point(antiIdealPoint));
		extrPoints.back()->setObjVector(k, extrPoints.back()->findMissingCoordinate(*h, k));
		extrPoints.back()->addActiveHyperplane(h);
		h->addVertex(); //extrPoints.back()
	}

	// create adjacency lists for extreme points
	std::list<Point*>::iterator it1;
	int k = 0;
	std::list<Point*>::iterator it2;
	int l;
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		l = 0;
		for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
			if (k != l) {
				(*it1)->addAdjacentPoint((*it2)->get_adress());
			}
			l++;
		}
		k++;
	}

	// add p artificial dominated plans (search cone)
	std::vector<double> nv(lp.get_p());
	for (int k = 1; k <= lp.get_p(); k++) {
		for (int l = 0; l < lp.get_p(); l++) {
			nv[l] = 0;
		}
		nv[k - 1] = 1;
		boundingBox[k - 1] = new Hyperplane(nv,antiIdealPoint[k - 1]);
	}

	int obj = 0;
	for (it1 = extrPoints.begin(); it1 != extrPoints.end(); ++it1) {
		for (int k = 1; k <= lp.get_p(); k++) {
			if (obj != k) {
				(*it1)->addActiveHyperplane(boundingBox[k - 1]);
			}
		}
		obj++;
	}
}

void LowerBoundSet::print() {

	std::cout << "\n============= Lower Bound Set ==============\n\n";

	std::cout << "Anti-ideal point: ( ";
	for (int k = 0; k < lp.get_p() - 1; k++) {
		std::cout << antiIdealPoint[k] << " , ";
	}
	std::cout << antiIdealPoint[lp.get_p() - 1] << " )\n";

	// Hyperplane representation
	std::cout << "\nHyperplane representation(" << facets.size() << " facets):\n";
	std::list<Hyperplane*>::iterator it;
	for (it = facets.begin(); it != facets.end(); ++it) {
		std::cout << "( ";
		for (int k = 0; k < lp.get_p()-1; k++) {
			std::cout << (*it)->get_normalVector(k) << " , ";
		}
		std::cout << (*it)->get_normalVector(lp.get_p()-1) << " ) " << " = " << (*it)->get_rhs() << "  -> defined by " << (*it)->get_nbDefiningPts() << " points" << std::endl;
	}

	// Extreme points representation

	std::list<Point*>::iterator it2;
	int s = 0;
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->isOnBoundingBox()) {
			++s;
		}
	}

	std::cout << "\n Extreme points representation (" << s << " extreme points + " << extrPoints.size() - s <<" artificial points):\n";
	for (it2 = extrPoints.begin(); it2 != extrPoints.end(); ++it2) {
		if (!(*it2)->isOnBoundingBox()) {
			std::cout << " ( ";
			for (int k = 0; k < lp.get_p() - 1; k++) {
				std::cout << (*it2)->get_objVector(k) << " , ";
			}
			std::cout << (*it2)->get_objVector(lp.get_p()-1) << " )\n";


			//int v = (int) sqrt(lp.get_n());
			//for (int i = 0; i < v; i++) {
			//	for (int j = 0; j < v; j++) {
			//		if ((*it2)->get_preImage(i * v + j) == 1) {
						//std::cout << lp.get;
			//		}
			//	}
			//	std::cout << "\n";
			//}
		}
	}
	std::cout << "\n";
}

void LowerBoundSet::compute(MathematicalFormulation& LP) {

	bool feasible(true);
	std::vector<double> y;
	Hyperplane* H;
	//int ctr = 0;

	// For exploring and updating list of extreme points at the same time
	std::list<Point*>::iterator checkpoint;
	std::list<Point*>::iterator currentPoint = extrPoints.begin();
	bool firstCheckpointReached = false;

	do
	{
		//print();
		//std::cout << "\n New point : ";
		//for (int i = 0; i < lp.get_p(); i++) {
		//	std::cout << " " << (*currentPoint)->get_objVector(i);
		//}
		nbIterations++;
		//print();
		//if (nbIterations >= 110){
		//	std::cout << "\n\n --------------- iteration " << nbIterations << " -----------------\n";
			//print();
		//}

		if (!warmstarted) {
			feasible = feasibilityModel->solveFeasibility(*(*currentPoint)->get_objVector());
		}
		else {
			// to do (case where LB is warmstarted) -> use point.feasible flag instead ???
		}
		//std::cout << " ---> feasible : " << feasible << std::endl;
		if (feasible) { // if this extreme point is a feasible objective vector, store preimage
			if (!warmstarted) { // preimage already computed if feasible when warmstarted
				(*currentPoint)->isNowFeasible(feasibilityModel);
			}					
			checkpoint = currentPoint;
			firstCheckpointReached = true;
			currentPoint++;
		}
		else // else, compute the cutting hyperplane
		{
			//toErase.push_back(bensonVertexIterator);
			y = bestValidPoint->getBestValidPoint(*(*currentPoint)->get_objVector(),interiorPoint);

			//std::cout << "\n\n best valid pts y = ";
			//for (int k = 0; k < LP.get_p(); k++) {
			//	std::cout << y[k] << " ";
			//}
			//std::cout << "\n";

			dualBenson->solveDualBenson(y);
			std::vector<double> normalVector = dualBenson->getNormalVector();
			double rhs = dualBenson->getRhs(LP);
			H = new Hyperplane(normalVector, rhs);
			//std::cout << "   -> Hyperplane found : ( ";
			//for (int k = 0; k < lp.get_p() - 1; k++) {
			//	std::cout << H->get_normalVector(k) << " , ";
			//}
			//std::cout << H->get_normalVector(lp.get_p() - 1) << " ) = " << rhs << std::endl;
			updatePolytope(H);

			// delete the non-feasible point
			//extrPoints.erase(currentPoint);
			if (firstCheckpointReached) {
				currentPoint = checkpoint;
				currentPoint++;
			}
			else
			{
				currentPoint = extrPoints.begin();
			}
		}

		// check if on bounding box

		// filter out-of-polytope hyperplanes using point ids

	} while (currentPoint != extrPoints.end());
	// --> include the last point to check completely ???? seems to be done actually
	// --> skip the first point ??? (anti-ideal)
	classifyExtremePoints();
}

void LowerBoundSet::updatePolytope(Hyperplane* H) {

	facets.push_back(H);
	//std::cout << "\n --- Update Polytope ---";

	//if (nbIterations >= 110) {
	//	std::cout << "Tested Hyperplane: ";
	//	H->print();
	//}

	// search for infeasible points
	std::list<Point*>::iterator it = extrPoints.begin();
	//std::list<Point*> discardedPts;
	while (it != extrPoints.end()) {
		(*it)->becomesNonDegenerate();
		if ((*it)->locatedOn(*H)) { // point below H, so it is not in the polytope anymore and it should be discarded   /// !H->above(**it)
			(*it)->becomesDegenerate();
			//discardedPts.push_back(it->get_adress());
			//if (nbIterations >= 110) {
			//	std::cout << "\npoint degenerated : ";
			//	(*it)->print();
			//}
		}
		else if ((*it)->below(*H)) {
			(*it)->becomesDiscarded();
			//if (nbIterations >= 110) {
			//	std::cout << "\npoint discarded : ";
			//	(*it)->print();
			//}
		}
		else {
			//std::cout << "\nnot this one...";
		}
		++it;
	}

	// search for new extreme points (comparison with edges, by adjacency lists)
	std::list<Point*>::iterator itAdj;
	std::list<Point*>* adjacentVertex;
	Point* newPts;
	std::vector<Point*> allNewVertices(0);

	for (it = extrPoints.begin(); it != extrPoints.end(); ++it) { // for each extreme point
		if ((*it)->isDiscarded()) { // if it is discarded

			//if (nbIterations >= 110) {
			//	std::cout << "\n -> new infeasible vertex :";
			//	(*it)->print();
			//}

			adjacentVertex = (*it)->get_adjList(); // we look at its adjacency list
			for (itAdj = adjacentVertex->begin(); itAdj != adjacentVertex->end(); ++itAdj) {
				
				//if (nbIterations >= 110) {
				//	std::cout << "\n       . test with :";
				//	(*itAdj)->print();
				//}

				if (!(*itAdj)->isDiscarded() && !(*itAdj)->isDegenerate()) { // for each non-discarded adjacent vertex, we compute the new point

					//if (nbIterations >= 110) {
					//	std::cout << "\n   -> Non-discarded adjecent vertex :";
					//	(*itAdj)->print();
					//}

					newPts = new Point((*itAdj)->edgeIntersection(**it, *H)); // new pts as intersection of H and edge
					(*itAdj)->replaceAdjVertex((*it)->get_adress(), newPts); // update adj vertex for feasible one
					newPts->addAdjacentPoint((*itAdj)->get_adress()); // init adj vertex for new one
					newPts->updateActiveHyperplanes(**it, **itAdj, H); // add a new active hyperplanes [ICI] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					extrPoints.push_back(newPts); //new pts added to list of extreme points
					allNewVertices.push_back(newPts); // remember the new vertex
					//H->addVertex(); //new point to hyperplane H

					//if (nbIterations >= 110) {
					//	std::cout << "\n      => new vertex : ";
					//	newPts->print();
					//}

					//std::cout << " with active hyperplanes -> ";
					//std::list<Hyperplane*>::iterator hp;
					//for (hp = newPts->get_activeHyperplanes()->begin(); hp != newPts->get_activeHyperplanes()->end(); hp++) {
					//	std::cout << *hp << " ";
					//}
				//}					
				}
			}
		}
		else if ((*it)->isDegenerate()) {
			(*it)->addActiveHyperplane(H);
			H->addVertex(); //new point to hyperplane H
			//if (abs(H->get_d() + 36.2143) <= 0.01) {
			//	std::cout << " Hyperplane : " << H;
			//	(*it)->print();
			//	std::cout << "\n";
			//}
			allNewVertices.push_back(*it); // remember the degenerate vertex
		}
	}

	// Update the adjacency lists for the new vertices + counters hyperplanes
	int s = allNewVertices.size();
	std::list <Hyperplane*>::iterator f;
	std::list<Hyperplane*>* hpp;
	//std::cout << "\n" << s << " new vertices !";
	//for (int i = 0; i < s; i++) {
	//	std::cout << "\n       new :";
	//	for (int k = 0; k < lp.get_p(); k++) {
	//		std::cout << " " << allNewVertices[i]->get_objVector(k);
	//	}
	//}
	
	for (int i = 0; i < s; i++) {
		for (int j = i + 1; j < s; j++) {
			//std::cout << "  ... we pass here!" << std::endl;
			if (allNewVertices[i]->sizeIntersectionActiveHyperplanes(*allNewVertices[j]) == H->get_dim()) { // is this true for p >= 3 ?
				allNewVertices[i]->addAdjacentPoint(allNewVertices[j]);
				allNewVertices[j]->addAdjacentPoint(allNewVertices[i]);
			}
		}
		//std::cout << "bwilibwi ?....";
		if (!allNewVertices[i]->isDegenerate()) {
			//std::cout << "bwalahaha\n";
			hpp = allNewVertices[i]->get_activeHyperplanes();
			for (f = hpp->begin(); f != hpp->end(); f++) {
				(*f)->addVertex();
				//if (abs((*f)->get_d() + 36.2143 ) <= 0.01) {
				//	std::cout << " Hyperplane : " << *f;
				//	allNewVertices[i]->print();
				//	std::cout << "\n";
				//}
			}
		}
	}

	// delete non-feasible vertices & update hyperplanes
	std::list<Point*>::iterator itdelete = extrPoints.begin();
	std::list<Point*>::iterator current;
	//std::list <Hyperplane*>::iterator f;
	std::list <Hyperplane*>::iterator currentface;
	std::list <Hyperplane*>* listHpp;

	while (itdelete != extrPoints.end()) {
		current = itdelete;
		++itdelete;
		if ((*current)->isDiscarded()) {
			listHpp = (*current)->get_activeHyperplanes();
			for (f = listHpp->begin(); f != listHpp->end(); f++) {
				(*f)->removeVertex();
			}
			(*current)->~Point(); // destroy this point
			extrPoints.erase(current);
		}
	}

	// delete redundants hyperplanes
	f = facets.begin();
	while (f != facets.end()) {
		currentface = f;
		++f;
		if ((*currentface)->isRedundant()) { //*H
			(*currentface)->~Hyperplane();
			facets.erase(currentface);
		}
	}
	//std::cout << "\n\n";
}

void LowerBoundSet::classifyExtremePoints() { // better name -> filterArtificialPoints ???

	for (std::list<Point*>::iterator pts = extrPoints.begin(); pts != extrPoints.end(); ++pts) {
		int k = 0;
		while (!(*pts)->isOnBoundingBox() && k < lp.get_p()) {
			if ((*pts)->get_objVector(k) == antiIdealPoint[k]) {
				(*pts)->becomesOnBoundingBox();
			}
			++k;
		}
	}

}

/* ==========================================================
		Getters
 ========================================================= */
